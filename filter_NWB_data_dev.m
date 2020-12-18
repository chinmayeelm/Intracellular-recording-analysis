filename = "M1_N3_T1";
filename_str = sprintf("%s.nwb", filename);
nwb_in = nwbRead(filename_str);

disp(nwb_in)


rec = nwb_in.acquisition.get('response_to_JO_stimulation');
data = rec.data.load;
rec_data = data(:,1);
stim_fb = data(:,2);
hes_data = data(:,3);
time = rec.timestamps.load;

no_of_protocols = 3;

no_of_trials = 5;
single_protocol_length = length(rec_data)/no_of_protocols;

single_trial_length = single_protocol_length/no_of_trials;

fs = 10000; %sampling freq

% time = 0:1/fs:24.0015;
% time = time(1:length(time)-1);



%% Run this if the order of the stimulus should be taken from a text file
[stim_order_sorted,idx] = sortfromtextfile(filename);


%% Run this to take stimulus order from nwb file stimulus description
[stim_order_sorted,idx] = sortfromnwb(nwb_in);

%%

rec_fc_low = 10;
rec_fc_high = 2000;
rec_filt_order = 8;
d_rec = designfilt('bandpassiir','FilterOrder',rec_filt_order, ...
'HalfPowerFrequency1',rec_fc_low,'HalfPowerFrequency2',rec_fc_high, ...
'SampleRate',fs, 'DesignMethod', 'butter');
% fvtool(d_rec);

hes_fc = 300;
[a_hes, b_hes] = butter(10, hes_fc/(fs/2));
 
    
filtered_data_bp = filtfilt(d_rec, rec_data);
hes_data_filtered = filtfilt(a_hes, b_hes, hes_data);
a= .9258; b=93.15; c=-1.455;
antennal_movement = (b./(hes_data_filtered -a)).^(1/3) + c;

fig_handle = consolidated_plot(time, filtered_data_bp, antennal_movement, stim_fb); 




%% Run this if the stimulus was not randomised
rec_protocols = reshape(filtered_data_bp, [no_of_protocols, single_protocol_length]);

rec_chir = filtered_data_bp(1:single_protocol_length);
rec_chirp = reshape(rec_chir, [no_of_trials, single_trial_length]);

rec_amp_sweep = rec_protocols(2,:);
rec_white_noise = rec_protocols(3,:);

%% Run this if the stimulus was randomised
rec_protocols_sorted = reshape_sort_data(filtered_data_bp, single_trial_length, no_of_protocols, no_of_trials, idx);
stim_protocols_hes_sorted = reshape_sort_data(antennal_movement, single_trial_length, no_of_protocols, no_of_trials, idx);
stim_protocols_ifb_sorted =  reshape_sort_data(stim_fb, single_trial_length, no_of_protocols, no_of_trials, idx);

%P = struct();
P = create_structs(time, rec_protocols_sorted,stim_protocols_hes_sorted,fs,stim_protocols_ifb_sorted,no_of_protocols, no_of_trials, single_trial_length, filename, stim_order_sorted);



%%

function P = create_structs(time, rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted, no_of_protocols, no_of_trials, single_trial_length, filename, stim_order_sorted)

    
    for i=1:no_of_protocols
       
        P(i).rec = rec_protocols_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        P(i).stim_hes = stim_protocols_hes_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        P(i).stim_ifb = stim_protocols_ifb_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        
        P(i).stim_name = stim_order_sorted(i*no_of_trials);
        
        [raster_data,gcfr]  = plotdata(no_of_trials, single_trial_length, fs, 100*P(i).rec, P(i).stim_hes, time, filename,  P(i).stim_name);
        P(i).raster = raster_data;
        P(i).gcfr = gcfr;

    end

end


function reshaped_sorted_data = reshape_sort_data(data, single_trial_length, no_of_protocols, no_of_trials, idx)
    reshaped_data = (reshape(data, [single_trial_length, no_of_protocols*no_of_trials]))';
    reshaped_sorted_data = reshaped_data(idx,:);
end


function [stim_order_sorted,idx] = sortfromnwb(nwb_in)
   
    stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
    stim_order = stim.stimulus_description;

    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);

end

function [stim_order_sorted,idx] = sortfromtextfile(filename)
    txt_file = sprintf("%s.txt", filename);
    fileID = fopen(txt_file, 'r');
    stim_order = fscanf(fileID, '%s');

    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);
end

function P = savefigures(filename, stim_name, figurehandle)
    
    png_name = sprintf("%s_%s.png", filename, stim_name);
    fig_name = sprintf("%s_%s.fig", filename, stim_name);
    
    savefig(figurehandle, fig_name);
    saveas(figurehandle, png_name);

end

function [raster_data,avg_gcfr]  = plotdata(no_of_trials, single_trial_length, fs, P_rec, P_stim_hes, time, filename, stim_name)
    
    

    raster_data = zeros(no_of_trials, single_trial_length);
    
    fig = figure();
    subplot(4,1,2);
    k=0.5;
    for i=1:no_of_trials
        [p,l] = findpeaks(P_rec(i,:), "MinPeakHeight",0.2*max(P_rec(i,:)));
        
        raster_data(i,l) = 1;
        spike_time = l/fs;
        for j = 1:length(spike_time)
            line([spike_time(j) spike_time(j)], [k k+0.5], 'Color', 'k');
        end
        %plot(spiketime_trials(i,:), '.');
        k = k+1;
    end
    
   % ylim([0 10]);
    A2 = gca;
    ylabel('Trials');
    
    [p,l] = findpeaks(P_rec(1,:), "MinPeakHeight",0.2*max(P_rec(i,:)));
    subplot(4,1,1); plot(time(1:single_trial_length), P_rec(1, :), l/fs, p, 'r.');
    A1 = gca;
    ylabel('Voltage (mV)');
    

    sum_of_spikes = sum(raster_data, 1);

    L = 5000;
    alpha = 8;
    gauss_win = gausswin(L, alpha);
%     gcfr = conv(sum_of_spikes, gauss_win);
%     gcfr = gcfr/no_of_trials;
    avg_gcfr = (filter(gauss_win, 1, sum_of_spikes))/no_of_trials;
    subplot(4,1,3); plot(time(1:single_trial_length), avg_gcfr, 'Color', [0.2,0.3,0.49]);
    A3 = gca;
    ylabel('Avg. GCFR');
    
    subplot(4,1,4); plot(time(1:single_trial_length), P_stim_hes(1, :), 'Color', [0.6, 0.2,0]);
    A4 = gca;
    ylabel('Indenter feedback voltage');
    %ylabel('Antennal movement');
    
    linkaxes([A1,A2,A3,A4], 'x');
    
    savefigures(filename, stim_name, fig);
end    


function fig_handle = consolidated_plot(time, filtered_data_bp, antennal_movement, stim_fb)
    fig_handle = figure();

    subplot(3,1,1);plot(time, filtered_data_bp*100, 'k');
    A1 = gca;
    title('Response to antennal movement')
    ylabel('Voltage (mV)')


    subplot(3,1,2); hold on; plot(time, antennal_movement, 'Color', [0.6, 0.2,0]);
    ylabel('Antennal movement (mm)')
    A2 = gca;

    subplot(3,1,3); hold on; plot(time, stim_fb,'Color', [0.2,0.3,0.49]);
    ylabel('Indenter feedback (V)')
    xlabel('time (s)')
    A3 = gca;

    linkaxes([A1 A2 A3], 'x');

    net_movement_antenna = max(antennal_movement) - min(antennal_movement)

end



