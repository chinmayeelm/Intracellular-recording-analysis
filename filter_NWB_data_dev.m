filename = "M1_N3_T1";
filename_str = sprintf("%s.nwb", filename);
nwb_in = nwbRead(filename_str);

disp(nwb_in);


rec = nwb_in.acquisition.get('response_to_JO_stimulation');
data = rec.data.load;
rec_data = data(:,1);
stim_fb = data(:,2);
hes_data = data(:,3);
time = rec.timestamps.load;
max_chirp_frq = 120;
amp_sweep_frq = 5;
blwgn_fc = 300;

fs = 10000; %sampling freq

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

no_of_trials = 5;
no_of_protocols = length(stim_order_vector)/no_of_trials;
ON_dur = 10;
OFF_dur =3;
single_protocol_length = length(rec_data)/no_of_protocols;
single_trial_length = single_protocol_length/no_of_trials;
%%

fig_handle = consolidated_plot(time, filtered_data_bp, antennal_movement, stim_fb); 

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
 
[gain,phase] = phase_plot(P, OFF_dur, ON_dur, fs);
for i=1:5 %no_of_protocols
%     P(i).crossings = zc(i,:);
        
      P(i).phase = phase(i);
      P(i).gain = gain(i);
end

%% Zero crossing and Tuning curve

for i = 1:no_of_protocols
        if P(i).stim_type == "frq"
%             zero_crossing(P(i).stim_ifb(1,:), P(i).norm_gcfr, fs, ON_dur, OFF_dur);
             zero_crossing(freq_chirp, P(i).norm_gcfr, fs, ON_dur, OFF_dur);
        end
end

%% STA of Band limited white Gaussian Noise

for i = 1:no_of_protocols
    if P(i).stim_type == "blwgn"
        [STA_freq, occurances] = STA_analysis(P(i).raster, P(i).antennal_movement, 0.1, fs);
    end
end
%%
function antennal_movement  = butter_filtfilt(data, fc, fs, order)
    
    [a, b] = butter(order, fc/(fs/2));

    filtered_data = filtfilt(a, b, data);
    a= .9258; b=93.15; c=-1.455;
    antennal_movement = (b./(filtered_data -a)).^(1/3) + c;
end

function stim_freq = fft_stim(stim_matrix, Fs, L)


%     [m,~]  = size(stim_matrix);
%     stim_freq = zeros(m, 1);
    
%     for i=1:m
%         stim_filt = sgolayfilt(stim_matrix(i,:), 3, 51);
        Y = fft(stim_matrix);

%         figure()
%         plot(stim_matrix(i, :));

        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
%         length(P1(2:end-1))

       
        F = Fs*(0:(L/2))/L;
        f = F(2:end-1);
%         length(f)
%         figure;
%         plot(f,P1(2:end-1)) 
%         title('Single-Sided Amplitude Spectrum of X(t)')
%         xlabel('f (Hz)')
%         ylabel('|P1(f)|')

        [~,loc] = max(P1(2:end-1));
        stim_freq = f(loc);
    
%     end

end

function P = create_structs(time, rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted, no_of_protocols, no_of_trials, single_trial_length, filename, stim_order_sorted, max_chirp_frq,amp_sweep_frq, blwgn_fc)

    
    for i=1:no_of_protocols
        
        
        P(i).rec = 100*( rec_protocols_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:));
        P(i).hes_data_unfilt = stim_protocols_hes_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        P(i).stim_ifb = stim_protocols_ifb_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        
        P(i).stim_name = stim_order_sorted(i*no_of_trials);
        type_frq = split(P(i).stim_name, '_');
        P(i).stim_type = type_frq(1);
        
        
        if (P(i).stim_type == "sin" ) || (P(i).stim_type == "sqr") || (P(i).stim_type == "step")
            P(i).stim_period = 1/str2num(type_frq(2));
            
        elseif P(i).stim_type == "frq"
            order=10;
            P(i).antennal_movement =  butter_filtfilt(P(i).hes_data_unfilt(1,:), 1.3*max_chirp_frq, fs, order);
            
        
            
        elseif P(i).stim_type == "blwgn"
            order=10;
            P(i).antennal_movement = butter_filtfilt(P(i).hes_data_unfilt(1,:), 1.3*blwgn_fc, fs, order);
            
        elseif P(i).stim_type == "amp"
            order=3;
            P(i).antennal_movement =  butter_filtfilt(P(i).hes_data_unfilt(1,:), 1.3*amp_sweep_frq, fs, order);
            
        end

        
        
        [raster,avg_gcfr,complete_trials]  = get_raster_gcfr(no_of_trials, P(i).rec, single_trial_length);
        P(i).raster = raster;
        P(i).avg_gcfr = avg_gcfr;
        P(i).complete_trials = complete_trials;
        
        P(i).norm_gcfr = [];
        P(i).crossings = [];

    end

end

function reshaped_data = reshape_data(data, single_trial_length, no_of_protocols, no_of_trials)
    reshaped_data = (reshape(data, [single_trial_length, no_of_protocols*no_of_trials]))';
    
end

function reshaped_sorted_data = sort_data(reshaped_data, idx)
    reshaped_sorted_data = reshaped_data(idx,:);
end    

function [stim_order_sorted,idx] = sortfromnwb(nwb_in)
   
    idx = 0;
    stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
    stim_order = stim.stimulus_description;

    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);

end

function [stim_order_vector, stim_order_sorted,idx] = sortfromtextfile(filename)
    txt_file = sprintf("%s.txt", filename);
    fileID = fopen(txt_file, 'r');
    stim_order = fscanf(fileID, '%s');

    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);
end

function savefigures(filename, stim_name, figurehandle)
    
    png_name = sprintf("%s_%s.png", filename, stim_name);
    fig_name = sprintf("%s_%s.fig", filename, stim_name);
    
    savefig(figurehandle, fig_name);
    saveas(figurehandle, png_name);

end

function [raster_data,avg_gcfr,no_of_true_trials]   = get_raster_gcfr(no_of_trials, P_rec, single_trial_length)
    

        no_of_true_trials = 0;
        raster_data = zeros(no_of_trials, single_trial_length);
        for i=1:no_of_trials
%             p=[]; l=[];
            [p,l] =  findpeaks(P_rec(i,:), "MinPeakHeight",0.2*max(P_rec(i,:)));

            if mode(p)<5
                continue;
            else
                raster_data(i,l) = 1;
                no_of_true_trials = no_of_true_trials+1;
            end
        end   
%         complete_trials = no_of_true_trials;
%         P.raster = raster_data;
        sum_of_spikes = sum(raster_data, 1);

        L = 1000;
        alpha = 4;
        gauss_win = gausswin(L, alpha);
        avg_gcfr = (filter(gauss_win, 1, sum_of_spikes))/no_of_true_trials;

end    

function P = plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P, gcfr_max)

    for i=1:no_of_protocols
        
%         fig = figure(i);
        [p,l] = findpeaks(P(i).rec(1,:), "MinPeakHeight",0.2*max(P(i).rec(1,:)));
%         A1 = subplot(4,1,1); plot(time(1:single_trial_length), P(i).rec(1, :), l/fs, p, 'r.');
%         ylabel('Voltage (mV)');
%         title((join(split(P(i).stim_name,"_")," ")) +" Hz");
        
%         A2 = subplot(4,1,2);
        k = 0.5;
        for j = 1:P(i).complete_trials
            l = find(P(i).raster(j,:)==1);
            spike_time = l/fs;
            for m = 1:length(spike_time)
%                 line([spike_time(m) spike_time(m)], [k k+0.5], 'Color', 'k');
            end
            k = k+1;
        end
%         ylabel('Trials');
            
        P(i).norm_gcfr = P(i).avg_gcfr/gcfr_max;
%         A3 = subplot(4,1,3); plot(time(1:single_trial_length), P(i).norm_gcfr, 'Color', [0.2,0.3,0.49]);
%         ylabel('Normalised GCFR');
        
%         A4 = subplot(4,1,4); plot(time(1:single_trial_length), P(i).antennal_movement(1, :), 'Color', [0.6, 0.2,0]);
   
    %     ylabel('Indenter feedback voltage');
%         ylabel('Antennal movement');
%         xlabel('time(s)');
        
%         linkaxes([A1,A2,A3,A4], 'x');
    
%         savefigures(filename, P(i).stim_name, fig);
    end
    
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

function [gain,phase_lag] = phase_plot(P, OFF_dur, ON_dur, fs)
   zc=0; 
   
  for i = 1:5
%     figure();
      %simulate a Schmidt trigger to ditect zero crossings. 

      ant_mov = P(i).antennal_movement(1,OFF_dur*fs:(OFF_dur+ON_dur)*fs);
%       mean_pos = mean(ant_mov);
      stim_filt = ant_mov; %sgolayfilt(ant_mov, 1, 5001);
      stim_filt = stim_filt-mean(stim_filt);


      
%       gcfr = sgolayfilt(P(i).norm_gcfr(OFF_dur*fs:(OFF_dur+ON_dur)*fs), 1,5001);
      gcfr = P(i).norm_gcfr(OFF_dur*fs:(OFF_dur+ON_dur)*fs);
      gcfr = gcfr-mean(gcfr);
%       A2 = subplot(3,1,2); plot(gcfr);
         
       [r,lags] = xcorr(gcfr,stim_filt, 'coeff');
%          figure();
%          A3 = subplot(3,1,3); stem(lags,r);
         
        [val, ind] = max(r);
        
        
        time_lag = abs(lags(ind))/fs;
        phase_lag(i) = (time_lag/P(i).stim_period)*360; %in degrees
        
        gain(i) = (max(stim_filt)-min(stim_filt))/(max(gcfr)-min(gcfr));
        
        figure();
        polarscatter(phase_lag, gain);
        title('Useless plot');

  end
  
end

function zero_crossing(stim, resp, fs, ON_dur, OFF_dur)

    stim = stim(OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    resp = resp(OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    mean_pos = mean(stim);

    locs = [];
    zc = [];
    for j = 2:length(stim)
      if (stim(j-1)<=mean_pos && stim(j)>=mean_pos && stim(j+1)>mean_pos)
          zc(j) = 1;
      end
    end
    [val,locs ]= find(zc==1); 

    A1 =  subplot(2,1,1); plot(stim);hold on; plot(locs, stim(locs), 'rx'); yline(mean_pos); hold off;
    
    stim_clips = [];
    resp_clips = [];
    stim_freq = 1 ./ (diff(locs) / fs);
    for k= 2:length(locs)
         stim_clips = stim(locs(k-1):locs(k)); %figure(); subplot(2,1,1); plot(stim_clips) ;
         resp_clips = resp(locs(k-1):locs(k));% subplot(2,1,2); plot(resp_clips);

         %stim_freq(k-1) = fft_stim(stim_clips, fs, length(stim_clips));
         max_FR(k-1) = max(resp_clips);

         [r,lags] = xcorr(stim_clips- mean(stim_clips), resp_clips-mean(resp_clips), 'coeff');
    end
    
    size(stim_freq)
    size(max_FR)
    A2 = subplot(2,1,2); plot(stim_freq, max_FR, 'o');
    
%     linkaxes([A1, A2],'x');
 
end

function [STA_freq, occurances] = STA_analysis(raster_data, stimulus, window, fs)

    [m,n] = size(raster_data);
    STA_freq = [];
    for i=1:m
        spike_locs = find(raster_data(i,:)==1);
        for j=1:length(spike_locs)
            if (spike_locs(j)-window*fs)< 0 
                continue;
            end
            STA = stimulus((spike_locs(j)-window*fs):spike_locs(j));
            STA_freq(i,j) = fft_stim(STA, fs, window*fs);
        end
        
    end
    histogram(STA_freq);
        

end







