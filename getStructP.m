function P = getStructP(dataDirectory,filename,clip_data_flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

expt_date = replace(dataDirectory, '-','.');
filenameParts = split(filename, '_');
mothId = filenameParts(1);
cd (join(['D:\Work\Recordings\' string(expt_date) '\' 'raw\' mothId '\'], ''))

filename_str = sprintf("%s.nwb", filename);
% filepath = join(['E:\Recordings\' string(expt_date) '\' 'raw\' mothId '\' filename_str], '');
nwb_in = nwbRead(filename_str); 

rec = nwb_in.acquisition.get('response_to_JO_stimulation');
data = rec.data.load;
time = rec.timestamps.load;
stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
intendedStimulus = stim.data.load;


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
gauss_win_sigma = 0.03; % 30 ms % fs is multiplied in the code later.

flag_meas_table = readtable('D:\Work\Recordings\Antenna flagellum measurements\flagellum-length-measurements.xlsx', 'VariableNamingRule','preserve');
flag_meas_table.Date = datetime(flag_meas_table.Date, 'format', 'dd-MM-uuuu');
flag_meas_table.Properties.RowNames = join([string(flag_meas_table.Date) flag_meas_table.MothID],'_');

filenameParts = split(filename,'_');
mothId = filenameParts(1);

rowId = join([string(datetime(replace(expt_date, '.','-'),'Format','dd-MM-uuuu')) mothId],'_');
try 
    movementRadius = table2array(flag_meas_table(rowId,'Lengthmm'));
    disp('flagellum measurement exists')

catch
    disp('Taking default radius')
    movementRadius = 0.64; %0.78; % in mm
end



stim_order_vector = string(split(stim.stimulus_description, ','));
[stim_order_sorted,idx] = sort(stim_order_vector);
no_of_protocols = length(unique(stim_order_sorted));
valid_trials = no_of_protocols * no_of_trials;


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
    time = time(1:length(data));
    valid_trials = round(length(data)/single_trial_length);
    
    stim_order_vector = stim_order_vector(1:valid_trials);
    [stim_order_sorted,idx] = sort(stim_order_vector);
    no_of_protocols = length(unique(stim_order_sorted));
    
    
elseif time(end)==0
    
    start_stim = OFF_dur*fs;
    stop_stim = (ON_dur+OFF_dur)*fs;
    single_trial_length = start_stim + stop_stim;
    if mod(length(data),single_trial_length)~=0
        single_trial_length = start_stim + stop_stim+1;
    end

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
    [stim_order_sorted,idx] = sort(stim_order_vector);
    no_of_protocols = length(unique(stim_order_sorted));

   
end
single_trial_length = length(data)/valid_trials;

rec_data = data(:,1);
stim_fb = data(:,2);
hes_data = data(:,3);

rec_fc_low = 10;
rec_fc_high = 2000;
rec_filt_order = 8; %8;
d_rec = designfilt('bandpassiir','FilterOrder',rec_filt_order, ...
    'HalfPowerFrequency1',rec_fc_low,'HalfPowerFrequency2',rec_fc_high, ...
    'SampleRate',fs, 'DesignMethod', 'butter');

filtered_data_bp = filtfilt(d_rec, rec_data);


% fig_handle = consolidated_plot(time, filtered_data_bp, hes_data, stim_fb, fs);
% title(join([expt_date replace(filename, '_', '')]));
 

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
d_ref = datetime('2021.06.01', 'InputFormat', 'yyyy.MM.dd');
expt_date = split(pwd, '\');
expt_date = expt_date(4);
d_check = datetime(expt_date, 'InputFormat', 'yyyy.MM.dd');

if d_check < d_ref
    a= .9258; b=93.15; c=-1.455; %before 1 June 2021
else
    a=0.5334; b=516.5; c = -3.233;
end

P = create_structs(rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted,intended_stimulus_sorted, no_of_protocols, no_of_trials, single_trial_length, stim_order_sorted,max_chirp_frq, amp_sweep_frq, blwgn_fc, ON_dur, a, b, c, gauss_win_L, gauss_win_sigma, movementRadius);

[P(:).ON_dur] = deal(ON_dur);
[P(:).OFF_dur] = deal(OFF_dur);
[P(:).time] = deal(time);
[P(:).filename] = deal(filename);
[P(:).date] = deal(expt_date);
[P(:).no_of_trials] = deal(no_of_trials);
[P(:).fs] = deal(fs);
[P(:).no_of_protocols] = deal(no_of_protocols);
[P(:).single_trial_length] = deal(single_trial_length);


end