filename = "M1_N2_T1";
filename_str = sprintf("%s.nwb", filename);
nwb_in = nwbRead(filename_str);

disp(nwb_in);


rec = nwb_in.acquisition.get('response_to_JO_stimulation');
data = rec.data.load;
rec_data = data(:,1);
stim_fb = data(:,2);
hes_data = data(:,3);
time = rec.timestamps.load;
% time = time(1:end-2);
max_chirp_frq =100;% 120;
amp_sweep_frq = 5;
blwgn_fc = 300;

ON_dur = 10;
OFF_dur =3;

no_of_trials = 5;

fs = 10000; %sampling freq

% stim_frequencies = find_stim_freq(stim_fb,ON_dur, OFF_dur, fs)

%% Run this if the order of the stimulus should be taken from a text file
[stim_order_vector, stim_order_sorted,idx] = sortfromtextfile(filename);


%% Run this to take stimulus order from nwb file stimulus description
% [stim_order_sorted,idx] = sortfromnwb(nwb_in);

%%

rec_fc_low = 10;
rec_fc_high = 2000;
rec_filt_order = 8;
d_rec = designfilt('bandpassiir','FilterOrder',rec_filt_order, ...
'HalfPowerFrequency1',rec_fc_low,'HalfPowerFrequency2',rec_fc_high, ...
'SampleRate',fs, 'DesignMethod', 'butter');
% fvtool(d_rec);

filtered_data_bp = filtfilt(d_rec, rec_data);

no_of_protocols = length(stim_order_vector)/no_of_trials;

single_protocol_length = length(rec_data)/no_of_protocols;
single_trial_length = single_protocol_length/no_of_trials;
%%

fig_handle = consolidated_plot(time, filtered_data_bp, hes_data, stim_fb); 

%% Run this if the stimulus was randomised
% [rec_protocols_reshaped, rec_protocols_sorted] = reshape_sort_data(filtered_data_bp, single_trial_length, no_of_protocols, no_of_trials, idx);
% [stim_protocols_hes_reshaped, stim_protocols_hes_sorted] = reshape_sort_data(antennal_movement, single_trial_length, no_of_protocols, no_of_trials, idx);
% [stim_protocols_ifb_reshaped,stim_protocols_ifb_sorted ]=  reshape_sort_data(stim_fb, single_trial_length, no_of_protocols, no_of_trials, idx);

rec_protocols_reshaped = reshape_data(filtered_data_bp, single_trial_length, no_of_protocols, no_of_trials);
stim_protocols_hes_reshaped = reshape_data(hes_data, single_trial_length, no_of_protocols, no_of_trials); %hes data not filtered. Antennal movement not calculated
stim_protocols_ifb_reshaped = reshape_data(stim_fb, single_trial_length, no_of_protocols, no_of_trials);

%% sort reshaped data

rec_protocols_sorted = sort_data(rec_protocols_reshaped,idx);
stim_protocols_hes_sorted = sort_data(stim_protocols_hes_reshaped, idx); %hes data not filtered. Antennal movement not calculated
stim_protocols_ifb_sorted = sort_data(stim_protocols_ifb_reshaped, idx);
%% Finding FFT - doesn't work if the duration of stimulus is less than the period

L = fs*ON_dur;
% 
stim_matrix = stim_protocols_hes_reshaped(:, fs*OFF_dur: fs*(OFF_dur+ON_dur));

stim_freq = fft_stim(stim_matrix, fs, L);

%% Create structs and plot data
P = create_structs(time, rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted, no_of_protocols, no_of_trials, single_trial_length, filename, stim_order_sorted,max_chirp_frq, amp_sweep_frq, blwgn_fc);

gcfr_all = extractfield(P, "avg_gcfr");
gcfr_max = max(gcfr_all);

%% Plot data

P = plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P, gcfr_max);

%% Phase plot

% Write code to select only sine, frq_chirp and amp_sweep protocols
% No :P


for i=1:no_of_protocols
%     P(i).crossings = zc(i,:);
    if P(i).stim_type == "sin"
        [gain,phase] = phase_plot(P(i).antennal_movement,P(i).norm_gcfr, OFF_dur, ON_dur, fs, P(i).stim_period);
        
      P(i).phase = phase;
      P(i).gain = gain;
    end
end

%% Tuning curve and spike phase Vs Frequency

for i = 1:no_of_protocols
        if P(i).stim_type == "frq"
%             zero_crossing(P(i).stim_ifb(1,:), P(i).norm_gcfr, fs, ON_dur, OFF_dur);
             [stim_freq, max_FR] = tuning_curve(P(i).antennal_movement(1,:), P(i).norm_gcfr, fs, ON_dur, OFF_dur); %require edit in the function to take all rows
             [I_spike_phase, II_spike_phase, III_spike_phase, stim_freq, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(P(i).antennal_movement(1,:), P(i).raster(1,:), fs, ON_dur, OFF_dur);
             figure(); scatter(I_spike_freq,I_spike_phase, 'rx'); hold on; scatter(II_spike_freq,II_spike_phase, 'bo');  scatter(III_spike_freq,III_spike_phase, 'k.')
             
             P(i).frq_chirp_f = stim_freq;
             P(i).frq_chirp_FR = max_FR;
             
             P(i).I_spike = [I_spike_freq', I_spike_phase'];
             P(i).II_spike = [II_spike_freq', II_spike_phase'];
             P(i).III_spike = [III_spike_freq', III_spike_phase'];
        end
end

%% STA of Band limited white Gaussian Noise

for i = 1:no_of_protocols
    if P(i).stim_type == "blwgn"
        [~, power_fft, frq_fft, STA]  = STA_analysis(P(i).raster, P(i).antennal_movement, 0.1, fs);
        P(i).power_fft = power_fft;
        P(i).frq_fft = frq_fft;
        P(i).STA = STA;
    end
end

%% Latency to spike in square or step wave stimulus

for i = 1:no_of_protocols
    
    if (P(i).stim_type == "sqr") || (P(i).stim_type == "step")
        [P(i).latencies, P(i).mean_latency] = spike_latency(P(i).antennal_movement,P(i).raster,ON_dur,OFF_dur,fs, P(i).complete_trials);
    end
        
end

mean_spike_latency = mean(extractfield(P, 'latencies'))
spike_latency_SD = std(extractfield(P, 'latencies'))
spike_latency_var = var(extractfield(P, 'latencies'))

%% Write structure to text file
analysis_data_file = sprintf("%s-analysis.txt", filename);
writetable(struct2table(P),analysis_data_file );


