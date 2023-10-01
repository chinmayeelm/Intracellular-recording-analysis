% clear;


expt_date = split(pwd, '\');
expt_date = expt_date(4);

expt_date = datetime(replace(expt_date, '.','-'),'Format','dd-MM-uuuu');


filename = "M1_N3_ramp";
filename_str = sprintf("%s.nwb", filename);
nwb_in = nwbRead(filename_str); 
clip_data_flag =0;

flag_meas_table = readtable('Antenna flagellum measurements\flagellum-length-measurements.xlsx');
flag_meas_table.Date = datetime(flag_meas_table.Date, 'format', 'dd-MM-uuuu');
filenameParts = split(filename,'_');
mothId = filenameParts(1);
flag_meas_table.Properties.RowNames = join([string(flag_meas_table.Date) flag_meas_table.MothID],'_');

rowId = join([string(expt_date) mothId],'_');
try 
    movementRadius = table2array(flag_meas_table(rowId,'Lengthmm'));
    disp('flagellum measurement exists')

catch
    disp('Taking default radius')
    movementRadius = 0.64; %0.78; % in mm
end

% disp(nwb_in);
rec = nwb_in.acquisition.get('response_to_JO_stimulation');
data = rec.data.load;
% rec_data = data(:,1);
% stim_fb = data(:,3); %2
% hes_data = data(:,2); %3

% time = linspace(0,length(data),1);
time = rec.timestamps.load;
stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
intendedStimulus = stim.data.load;

%sampling freq

max_chirp_frq =  150;
max_sqr_sin_frq = 10;
amp_sweep_frq = 5;
blwgn_fc = 300;

parameters = nwb_in.general_stimulus.load;
if length(parameters) == 5
   
    ON_dur = str2num(parameters(3));
    OFF_dur =  str2num(parameters(4));
    no_of_trials = str2num(parameters(5));
    fs = str2num(parameters(1));

elseif length(parameters) == 6
    
    ON_dur = str2num(parameters(4));
    OFF_dur =  str2num(parameters(5));
    no_of_trials = str2num(parameters(6));
    fs = str2num(parameters(1));
    
else
    
    disp('No parameters in the file');
end

gauss_win_L = fs/5;
gauss_win_sigma = 0.03; % 30 ms



% ON_dur = 4;
% OFF_dur =3;
% no_of_trials = 5;
% fs = 1e4;




% stim_frequencies = find_stim_freq(stim_fb,ON_dur, OFF_dur, fs)

% Run this if the order of the stimulus should be taken from a text file
% [stim_order_vector, stim_order_sorted,idx] = sortfromtextfile(filename);


% Run this to take stimulus order from nwb file stimulus description
stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
%         stim_order = stim.stimulus_description;
stim_order_vector = string(split(stim.stimulus_description, ','));
% [stim_order_sorted, idx] = sortfromnwb(stim_order_vector);
% stim_order_vector = LUT_intra.stim_order(2);
[stim_order_sorted,idx] = sort(stim_order_vector);
no_of_protocols = length(unique(stim_order_sorted));
valid_trials = no_of_protocols * no_of_trials;

% clear nwb_in

start_stim = OFF_dur*fs;
stop_stim = (ON_dur+OFF_dur)*fs;

if clip_data_flag == 1
    
    start_stim = OFF_dur*fs;
    stop_stim = (ON_dur+OFF_dur)*fs;
    single_trial_length = start_stim + stop_stim;
    if mod(length(data),single_trial_length)~=0
        single_trial_length = start_stim + stop_stim+1;
    end

    prompt1 = 'Enter start point for data';
    prompt2 = 'Enter stop point for data';
    
    start_point = input(prompt1)*fs;
    stop_point = input(prompt2)*fs;
    start_clip_point = start_point-mod(start_point,single_trial_length)+1;
    stop_clip_point = stop_point-mod(stop_point,single_trial_length);
    
    
    data = data(start_clip_point:stop_clip_point,:);
%     time = time(start_clip_point:stop_clip_point);
    time = time(1:length(data));
    valid_trials = round(length(data)/single_trial_length);
    
    start_valid_trial = ((start_clip_point-1)/single_trial_length)+1;
    stop_valid_trial = stop_clip_point/single_trial_length;
    stim_order_vector = stim_order_vector(start_valid_trial:stop_valid_trial);
    intendedStimulus = intendedStimulus(:,start_valid_trial:stop_valid_trial); 
    [stim_order_sorted,idx] = sort(stim_order_vector);
    no_of_protocols = length(unique(stim_order_sorted));
    
elseif time(end)==0
    
    start_stim = OFF_dur*fs;
    stop_stim = (ON_dur+OFF_dur)*fs;
    single_trial_length = start_stim + stop_stim;
    if mod(length(data),single_trial_length)~=0
        single_trial_length = start_stim + stop_stim+1;
    end

    stop_point = find(data(:,2)==0,1) - 1;
    clip_point = stop_point-mod(stop_point,single_trial_length);
    
    data = data(1:clip_point,:);
    time = time(1:clip_point);
    valid_trials = round(length(data)/single_trial_length);
    
    stim_order_vector = stim_order_vector(1:valid_trials);
    intendedStimulus = intendedStimulus(:,1:valid_trials); 
    [stim_order_sorted,idx] = sort(stim_order_vector);
    no_of_protocols = length(unique(stim_order_sorted));
end
%
single_trial_length = length(data)/valid_trials;
%     no_of_protocols = length(stim_order_sorted)/no_of_trials;

% valid_trials = round(length(data)/single_trial_length);

rec_data = data(:,1);
stim_fb = data(:,2);
hes_data = data(:,3);

rec_fc_low = 10;
rec_fc_high = 2000;
rec_filt_order = 8; %8;
d_rec = designfilt('bandpassiir','FilterOrder',rec_filt_order, ...
    'HalfPowerFrequency1',rec_fc_low,'HalfPowerFrequency2',rec_fc_high, ...
    'SampleRate',fs, 'DesignMethod', 'butter');
% fvtool(d_rec);

filtered_data_bp = filtfilt(d_rec, rec_data);

% no_of_protocols = length(stim_order_sorted)/no_of_trials;

% single_protocol_length = length(rec_data)/no_of_protocols;
% single_trial_length = single_protocol_length/no_of_trials;
%

fig_handle = consolidated_plot(time, filtered_data_bp, hes_data, stim_fb, fs);
% 

% Run this if the stimulus was randomised
% [rec_protocols_reshaped, rec_protocols_sorted] = reshape_sort_data(filtered_data_bp, single_trial_length, no_of_protocols, no_of_trials, idx);
% [stim_protocols_hes_reshaped, stim_protocols_hes_sorted] = reshape_sort_data(antennal_movement, single_trial_length, no_of_protocols, no_of_trials, idx);
% [stim_protocols_ifb_reshaped,stim_protocols_ifb_sorted ]=  reshape_sort_data(stim_fb, single_trial_length, no_of_protocols, no_of_trials, idx);

rec_protocols_reshaped = reshape_data(filtered_data_bp, single_trial_length, no_of_protocols, valid_trials);
stim_protocols_hes_reshaped = reshape_data(hes_data, single_trial_length, no_of_protocols, valid_trials); %hes data not filtered. Antennal movement not calculated
stim_protocols_ifb_reshaped = reshape_data(stim_fb, single_trial_length, no_of_protocols, valid_trials);

% sort reshaped data

rec_protocols_sorted = sort_data(rec_protocols_reshaped,idx);
stim_protocols_hes_sorted = sort_data(stim_protocols_hes_reshaped, idx); %hes data not filtered. Antennal movement not calculated
stim_protocols_ifb_sorted = sort_data(stim_protocols_ifb_reshaped, idx);
intended_stimulus_sorted = sort_data(intendedStimulus', idx);
% Clip 2s of baseline activity in the beginning and end of trials
%{
start_clip_point = 2*fs+1;
stop_clip_point = single_trial_length - 2*fs;

rec_protocols_sorted = rec_protocols_sorted(:,start_clip_point : stop_clip_point);
stim_protocols_hes_sorted = stim_protocols_hes_sorted(:,start_clip_point : stop_clip_point); %hes data not filtered. Antennal movement not calculated
stim_protocols_ifb_sorted = stim_protocols_ifb_sorted(:,start_clip_point : stop_clip_point);

single_trial_length = length(rec_protocols_sorted);
%}

% Create structs
d_ref = datetime('2021.06.01', 'InputFormat', 'yyyy.MM.dd') ;
d_check = datetime(expt_date, 'InputFormat', 'yyyy.MM.dd');

if d_check < d_ref
    a= .9258; b=93.15; c=-1.455; %before 1 June 2021
else
    a=0.5334; b=516.5; c = -3.233;
end

P = create_structs(rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted,intended_stimulus_sorted, no_of_protocols, no_of_trials, single_trial_length, stim_order_sorted,max_chirp_frq, amp_sweep_frq, blwgn_fc, ON_dur, a, b, c, gauss_win_L, gauss_win_sigma, movementRadius);
% end
%

%     gcfr_all = extractfield(P,"gcfr");
%     gcfr_max = max(gcfr_all);

[P(:).ON_dur] = deal(ON_dur);
[P(:).OFF_dur] = deal(OFF_dur);
[P(:).time] = deal(time);
[P(:).filename] = deal(filename);
[P(:).date] = deal(expt_date);
[P(:).no_of_trials] = deal(no_of_trials);
[P(:).fs] = deal(fs);
[P(:).no_of_protocols] = deal(no_of_protocols);
[P(:).single_trial_length] = deal(single_trial_length);

for i=1:no_of_protocols
    gcfr_max = max(P(i).gcfr(:));
    P(i).norm_gcfr = P(i).gcfr/gcfr_max;
    
end

% Plot data
% figure;
% plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P);

%%


% Phase plot             Not working

% Write code to select only sine, frq_chirp and amp_sweep protocols
% No :P

%{
for i=1:no_of_protocols
%     P(i).crossings = zc(i,:);
    if P(i).stim_type == "sin" || P(i).stim_type == "amp"
        [gain,phase] = phase_plot(P(i).antennal_movement,P(i).norm_gcfr, OFF_dur, ON_dur, fs, P(i).stim_period);
        
      P(i).phase = phase;
      P(i).gain = gain;
    end
end
%}
% GCFR Vs frequency and spike phase Vs Frequency

for i = 1:no_of_protocols
            if P(i).stim_type == "frq" || P(i).stim_type =="dec"
    
    
                 [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(P(i).antennal_movement(1,:), P(i).raster(1,:), fs, ON_dur, OFF_dur);
%                  figure(); scatter(I_spike_freq,I_spike_phase,100, 'k.'); hold on;
%                  pmin = min(I_spike_phase);
%                  pmax = max(I_spike_phase);
%                  pimin = floor(pmin/pi);
%                  pimax = ceil(pmax/pi);
%     %              yticks(0:pi/4:2*pi);
%                  yticks((pimin:pimax) * pi);
%                  yticklabels( string(pimin:pimax) + "\pi" )
%                  scatter(II_spike_freq,II_spike_phase, 100,  'b.');
%                  scatter(III_spike_freq,III_spike_phase, 100, 'r.')
%                  ylabel('Spike phase (rad)');
%                  xlabel('Stimulus frequency (Hz)');
%                  legend('I spike','II spike', 'III spike');
            end
    
    if P(i).stim_type == "frq"
        
        %              P(i).I_spike = [I_spike_freq', I_spike_phase'];
        %              P(i).II_spike = [II_spike_freq', II_spike_phase'];
        %              P(i).III_spike = [III_spike_freq', III_spike_phase'];
        
        inc_frq_chirp_f = linspace(1,max_chirp_frq,ON_dur*fs+1);
        P(i).inc_frq_chirp_f = inc_frq_chirp_f;
        figure;
        [lineOut, ~] = stdshade(P(i).gcfr(:,start_stim:stop_stim),0.2,'k',P(i).inc_frq_chirp_f);
        inc_chirp_gcfr = P(i).gcfr(:,start_stim:stop_stim);
        lineOut.LineWidth  = 0.05;
        lineOut.LineWidth  = 0.01;
        ylabel 'Firing rate (Hz)';
        xlabel 'Frequency (Hz)';
        title ('Response to increasing frequency chirp');

    elseif P(i).stim_type == "dec"
        P(i).dec_frq_chirp_f = linspace(max_chirp_frq,1,ON_dur*fs+1);
        
%         figure;
%         [lineOut, ~] = stdshade(P(i).gcfr(:,start_stim:stop_stim),0.2,'k',P(i).dec_frq_chirp_f);
%         lineOut.LineWidth  = 0.05;
%         ylabel 'Firing rate (Hz)';
%         xlabel 'Frequency (Hz)';
%         title ('Response to decreasing frequency chirp');
        %
        %              P(i).dec_chirp_I_spike = [I_spike_freq', I_spike_phase'];
        %              P(i).dec_chirp_II_spike = [II_spike_freq', II_spike_phase'];
        %              P(i).dec_chirp_III_spike = [III_spike_freq', III_spike_phase'];
    end
end

% STA of Band limited white Gaussian Noise

% for i = 1:no_of_protocols
%     if P(i).stim_type == "blwgn"
%         STA_window = 0.1;
%         STA  = STA_analysis(P(i).raster, P(i).antennal_movement, STA_window, fs, start_stim, stop_stim);
% 
%         [stim_freq, power_fft, frq_fft] = get_fft(STA, fs, window*fs);        
%         P(i).power_fft = power_fft;
%         P(i).frq_fft = frq_fft;
%         P(i).STA = STA;
%     end
% end



% Covariance analysis

for i = 1:no_of_protocols
    if P(i).stim_type == "blwgn"
        window = 0.03;
        [ev1, ev2, STA] = cov_analysis(P(i).raster, P(i).antennal_movement, window, fs, start_stim, stop_stim);
%         P(i).cov_matrix = cov_matrix;
                P(i).ev1 = ev1;
                P(i).ev2 = ev2;
                P(i).STA = STA;
        [power_fft, frq_fft] = get_fft(STA, fs, window*fs);        
        P(i).power_fft = power_fft;
        P(i).frq_fft = frq_fft;
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

% meta_struct(struct_num).moth = P;
% amp_sweep(struct_num) = P(1);
% blwgn(struct_num) = P(2);
% dec_frq_chirp(struct_num) = P;
% inc_frq_chirp(struct_num) = P;
% T_sqr_27_01(struct_num,:) = struct2table(P, 'AsArray', true);
% struct_num = struct_num+1;


