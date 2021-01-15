[~,yymm,dd] = fileparts(pwd);
date = strcat(yymm, dd);

filename = "M1_N1_T2";
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
% rec_data = data(1:end-2,1);
% stim_fb = data(1:end-2,2);
% hes_data = data(1:end-2,3);

max_chirp_frq =  120;
max_sqr_sin_frq = 120;
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


%% Create structs 
P = create_structs(rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted, no_of_protocols, no_of_trials, single_trial_length, stim_order_sorted,max_chirp_frq, amp_sweep_frq, blwgn_fc, max_sqr_sin_frq);

% gcfr_all = extractfield(P, "avg_gcfr");
gcfr_all = extractfield(P, "gcfr");
gcfr_max = max(gcfr_all);

[P(:).ON_dur] = deal(ON_dur);
[P(:).OFF_dur] = deal(OFF_dur);
[P(:).time] = deal(time);
[P(:).filename] = deal(filename);
[P(:).date] = deal(date);
[P(:).no_of_trials] = deal(no_of_trials);
[P(:).fs] = deal(fs);
[P(:).no_of_protocols] = deal(no_of_protocols);
[P(:).single_trial_length] = deal(single_trial_length);

for i=1:no_of_protocols
%     P(i).norm_gcfr = P(i).avg_gcfr/gcfr_max;
     P(i).norm_gcfr = P(i).gcfr/gcfr_max;

end

%% Plot data

P = plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P);

%% Phase plot             Not working

% Write code to select only sine, frq_chirp and amp_sweep protocols
% No :P


for i=1:no_of_protocols
%     P(i).crossings = zc(i,:);
    if P(i).stim_type == "sin" || P(i).stim_type == "amp"
        [gain,phase] = phase_plot(P(i).antennal_movement,P(i).norm_gcfr, OFF_dur, ON_dur, fs, P(i).stim_period);
        
      P(i).phase = phase;
      P(i).gain = gain;
    end
end

%% Tuning curve and spike phase Vs Frequency

for i = 1:no_of_protocols
        if P(i).stim_type == "frq"
%             zero_crossing(P(i).stim_ifb(1,:), P(i).norm_gcfr, fs, ON_dur, OFF_dur);
             [stim_freq, max_FR] = tuning_curve(P(i).antennal_movement, P(i).norm_gcfr, fs, ON_dur, OFF_dur); %require edit in the function to take all rows
             [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(P(i).antennal_movement(1,:), P(i).raster(1,:), fs, ON_dur, OFF_dur);
             figure(); scatter(I_spike_freq,I_spike_phase,100, 'k.'); hold on; 
             scatter(II_spike_freq,II_spike_phase, 100,  'b.');  
             scatter(III_spike_freq,III_spike_phase, 100, 'r.')
             ylabel('Spike phase (rad)');
             xlabel('Stimulus frequency (Hz)');
             legend('I spike','II spike', 'III spike');
             
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
        STA_window = 0.1;
        [~, power_fft, frq_fft, STA]  = STA_analysis(P(i).raster, P(i).antennal_movement, STA_window, fs);
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



%% Write structure to another structure

meta_struct(struct_num).moth = P;
struct_num = struct_num+1;


