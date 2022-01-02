% [~,yymm,dd] = fileparts(pwd);
% expt_date = strcat(yymm, dd);
% meta_table = [];
% LUT_intra.expt_date = datetime(LUT_intra.expt_date, 'Format', 'yyyy.MM.dd');

% chirp_meta_num = 1;
% ampsweep_meta_num = 1;
% blwgn_meta_num = 1;
% sin_meta_num = 1;
% sqr_meta_num = 1;
% impulse_meta_num = 1;

clip_data_flag = 0;
for LUT_intra_row_idx = 62
    
     
    
    
    LUT_intra_row_idx
    cd ..
    cd(string(LUT_intra.expt_date(LUT_intra_row_idx)))
    expt_date = LUT_intra.expt_date(LUT_intra_row_idx);
    filename = LUT_intra.filename(LUT_intra_row_idx);
    
%     filename = "M3_N3_T3"; %change HES parameters if required.
    filename_str = sprintf("%s.nwb", filename);
    nwb_in = nwbRead(filename_str);
    
    % disp(nwb_in);
    rec = nwb_in.acquisition.get('response_to_JO_stimulation');
    data = rec.data.load;
    time = rec.timestamps.load;
    
%     data = downsample(data,10);
%     time = downsample(time,10);
    
    %temporary parameters to check files not on LUT
%     parameters = nwb_in.general_stimulus.load;
%     ON_dur = str2num(parameters(3));
%     OFF_dur =str2num(parameters(4));
%     no_of_trials = str2num(parameters(5));
%     fs = str2num(parameters(1));
%     max_chirp_frq = 150;
%     amp_sweep_frq = 5;
%     blwgn_fc = 300; 
    
    
    row_name = strcat(string(expt_date), "_",filename);
    
    parameters = LUT_intra(row_name,:);     
    ON_dur = parameters.ON_dur;
    OFF_dur = parameters.OFF_dur;
    no_of_trials = parameters.no_of_trials;
    fs = parameters.fs;
    max_chirp_frq = parameters.max_chirp_frq;
    amp_sweep_frq = parameters.amp_sweep_frq;
    blwgn_fc = parameters.blwgn_fc;
    impulse_dur = ON_dur;
    
    
    start_stim = OFF_dur*fs;
    stop_stim = (ON_dur+OFF_dur)*fs;
    
    single_trial_length = start_stim + stop_stim+1;
    
    
    % Run this if the order of the stimulus should be taken from a text file
    % [stim_order_vector, stim_order_sorted,idx] = sortfromtextfile(filename);
    
    
    % Run this to take stimulus order from nwb file stimulus description
    
    if  parameters.stim_order == ""
        stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
        %         stim_order = stim.stimulus_description;
        stim_order_vector = string(split(stim.stimulus_description, ','));
    else
        stim_order = parameters.stim_order;
        stim_order_vector = split(parameters.stim_order, ',');
    end
    
    
    
    if time(end)==0
        
        stop_point = find(data(:,2)==0,1) - 1;
        clip_point = stop_point-mod(stop_point,single_trial_length);
        
        data = data(1:clip_point,:);
        time = time(1:clip_point);
        valid_trials = length(data)/single_trial_length;
        
        stim_order_vector = stim_order_vector(1:valid_trials);
        
    elseif clip_data_flag == 1
        
        prompt = 'Enter stop point for data';
        stop_point = input(prompt)*fs;
        clip_point = stop_point-mod(stop_point,single_trial_length);
        
        data = data(1:clip_point,:);
        time = time(1:clip_point);
        valid_trials = length(data)/single_trial_length;
        
        stim_order_vector = stim_order_vector(1:valid_trials);
    end
    
    %     [stim_order_sorted,stim_order_vector, idx] = sortfromnwb(stim_order);
    [stim_order_sorted,idx] = sort(stim_order_vector);
    %     no_of_protocols = length(stim_order_sorted)/no_of_trials;
    no_of_protocols = length(unique(stim_order_sorted));
    valid_trials = round(length(data)/single_trial_length);
    
    %     single_protocol_length = length(rec_data)/no_of_protocols;
    %     single_trial_length = single_protocol_length/no_of_trials;
    
    rec_data = data(:,1);
    stim_fb = data(:,2);
    hes_data = data(:,3);
    
    rec_fc_low = 10;
    rec_fc_high = 2000;
    rec_filt_order = 8;
    d_rec = designfilt('bandpassiir','FilterOrder',rec_filt_order, ...
        'HalfPowerFrequency1',rec_fc_low,'HalfPowerFrequency2',rec_fc_high, ...
        'SampleRate',fs, 'DesignMethod', 'butter');
    % fvtool(d_rec);
    
    filtered_data_bp = filtfilt(d_rec, rec_data);
    
    
    %
    
%     fig_handle = consolidated_plot(time, filtered_data_bp, hes_data, stim_fb, fs);
    
    %     pause;
    %
    % Run this if the stimulus was randomised
    % [rec_protocols_reshaped, rec_protocols_sorted] = reshape_sort_data(filtered_data_bp, single_trial_length, no_of_protocols, no_of_trials, idx);
    % [stim_protocols_hes_reshaped, stim_protocols_hes_sorted] = reshape_sort_data(antennal_movement, single_trial_length, no_of_protocols, no_of_trials, idx);
    % [stim_protocols_ifb_reshaped,stim_protocols_ifb_sorted ]=  reshape_sort_data(stim_fb, single_trial_length, no_of_protocols, no_of_trials, idx);
    
    %     rec_protocols_reshaped = reshape_data(filtered_data_bp, single_trial_length, no_of_protocols, no_of_trials);
    %     stim_protocols_hes_reshaped = reshape_data(hes_data, single_trial_length, no_of_protocols, no_of_trials); %hes data not filtered. Antennal movement not calculated
    %     stim_protocols_ifb_reshaped = reshape_data(stim_fb, single_trial_length, no_of_protocols, no_of_trials);
    
    rec_protocols_reshaped = reshape_data(filtered_data_bp, single_trial_length, no_of_protocols, valid_trials);
    stim_protocols_hes_reshaped = reshape_data(hes_data, single_trial_length, no_of_protocols, valid_trials); %hes data not filtered. Antennal movement not calculated
    stim_protocols_ifb_reshaped = reshape_data(stim_fb, single_trial_length, no_of_protocols, valid_trials);
    
    % sort reshaped data
    
    %     rec_protocols_sorted = sort_data(rec_protocols_reshaped,idx);
    %     stim_protocols_hes_sorted = sort_data(stim_protocols_hes_reshaped, idx); %hes data not filtered. Antennal movement not calculated
    %     stim_protocols_ifb_sorted = sort_data(stim_protocols_ifb_reshaped, idx);
    rec_protocols_sorted = rec_protocols_reshaped(idx,:);
    stim_protocols_hes_sorted = stim_protocols_hes_reshaped(idx,:); %hes data not filtered. Antennal movement not calculated
    stim_protocols_ifb_sorted = stim_protocols_ifb_reshaped(idx,:);
    
    % Create structs
    
    d_ref = datetime('2021.06.01', 'InputFormat', 'yyyy.MM.dd');
    expt_date = split(pwd, '\');
    expt_date = expt_date(3);
    d_check = datetime(expt_date, 'InputFormat', 'yyyy.MM.dd');
    
    if d_check < d_ref
        a= .9258; b=93.15; c=-1.455; %before 1 June 2021
    else
        a=0.5334; b=516.5; c = -3.233;
    end
    
    P = create_structs(rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted, no_of_protocols, no_of_trials, single_trial_length, stim_order_sorted,max_chirp_frq, amp_sweep_frq, blwgn_fc, ON_dur, a, b, c);
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
    
%     plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P);

%     pause;
%    
    % Phase plot
end

%%
    
    for i=2:no_of_protocols
        %     P(i).crossings = zc(i,:);
        if P(i).stim_type == "sin" %|| P(i).stim_type == "amp"
            [corr_coeff, vs, phase] = phase_plot(P(i).antennal_movement,P(i).norm_gcfr, P(i).raster, start_stim, stop_stim, fs, P(i).stim_period);
            
            P(i).phase = phase;
            P(i).vs = vs;
        end
    end
    
    % Amplitude Vs spike phase
    for i = 1:no_of_protocols
        if P(i).stim_type == "amp"
            
            
            [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_amp, II_spike_amp, III_spike_amp] = spike_phase(P(i).antennal_movement(1,:), P(i).raster(1,:), fs, start_stim, stop_stim, P(i).stim_type, 1/P(i).stim_period);
            figure(); scatter(I_spike_amp,I_spike_phase,100, '.k'); hold on;
            pmin = min([I_spike_phase II_spike_phase III_spike_phase]);
            pmax = max([I_spike_phase II_spike_phase III_spike_phase]);
            pimin = floor(pmin/pi);
            pimax = ceil(pmax/pi);
            ylim([pimin pimax]);
            yticks(pimin:pimax/4:pimax)
            yticklabels( string(pimin:pimax/4:pimax) + "\pi" )
            scatter(II_spike_amp,II_spike_phase, 100,  'b.');
            scatter(III_spike_amp,III_spike_phase, 100, 'r.');
            title('Amplitude sweep');
            ylabel('Spike phase (rad)');
            xlabel('Stimulus amplitude (mm)');
            box off;
            legend('I spike','II spike', 'III spike', 'Location', 'best');
            
            P(i).I_spike = [I_spike_amp', I_spike_phase'];
            P(i).II_spike = [II_spike_amp', II_spike_phase'];
            P(i).III_spike = [III_spike_amp', III_spike_phase'];
        end
    end
    % GCFR Vs frequency and spike phase Vs Frequency
    
    for i = 1:no_of_protocols
        if P(i).stim_type == "frq"  %|| P(i).stim_type =="dec"
            
            P(i).inc_frq_chirp_f = linspace(1,max_chirp_frq,ON_dur*fs+1);
%             [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(P(i).antennal_movement(1,:), P(i).raster(1,:), fs, start_stim, stop_stim, P(i).stim_type, P(i).inc_frq_chirp_f);
            [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(P(i).stim_ifb(1,:), P(i).raster(1,:), fs, start_stim, stop_stim, P(i).stim_type, P(i).inc_frq_chirp_f);
            figure(); scatter(I_spike_freq,I_spike_phase,100, '.k'); hold on;
            pmin = min([I_spike_phase II_spike_phase III_spike_phase]);
            pmax = max([I_spike_phase II_spike_phase III_spike_phase]);
            pimin = floor(pmin/pi);
            pimax = ceil(pmax/pi);
            ylim([pimin pimax]);
            yticks(pimin:pimax/4:pimax)
            yticklabels( string(pimin:pimax/4:pimax) + "\pi" )
            scatter(II_spike_freq,II_spike_phase, 100,  'b.');
            scatter(III_spike_freq,III_spike_phase, 100, 'r.');
            title('Frequency chirp');
            ylabel('Spike phase (rad)');
            xlabel('Stimulus frequency (Hz)');
            box off;
            legend('I spike','II spike', 'III spike', 'Location', 'best');
        end
        
        if P(i).stim_type == "frq"
            
            P(i).I_spike = [I_spike_freq', I_spike_phase'];
            P(i).II_spike = [II_spike_freq', II_spike_phase'];
            P(i).III_spike = [III_spike_freq', III_spike_phase'];
            
            P(i).inc_frq_chirp_f = linspace(0,max_chirp_frq,ON_dur*fs+1); % f0 = 0 for data 1 Sep onwards
            figure;
            [lineOut, ~] = stdshade(P(i).norm_gcfr(:,start_stim:stop_stim),0.2,'k',P(i).inc_frq_chirp_f);
            lineOut.LineWidth  = 0.01;
            ylabel 'GCFR';
            xlabel 'Frequency (Hz)';
            title ('Response to increasing frequency chirp');
            box off;
            
        elseif P(i).stim_type == "dec"
            P(i).dec_frq_chirp_f = linspace(max_chirp_frq,0,ON_dur*fs+1); % f0 = 0 for data 1 Sep onwards
            
            %             hold on;
            figure;
            [lineOut, ~] = stdshade(P(i).gcfr(:,start_stim:stop_stim),0.2,'k',fliplr(P(i).dec_frq_chirp_f));
            lineOut.LineWidth  = 0.05;
            ax = gca;
            ax.XTickLabel = flipud(ax.XTickLabel);
            ylabel 'GCFR';
            xlabel 'Frequency (Hz)';
            box off;
            title ('Response to decreasing frequency chirp');
            
            P(i).dec_chirp_I_spike = [I_spike_freq', I_spike_phase'];
            P(i).dec_chirp_II_spike = [II_spike_freq', II_spike_phase'];
            P(i).dec_chirp_III_spike = [III_spike_freq', III_spike_phase'];
        end
    end
    
    %     pause;
    
%%  
    % STA of Band limited white Gaussian Noise
    
    for i = 1:no_of_protocols
        if P(i).stim_type == "blwgn"% || P(i).stim_type == "frq"
            STA_window = 0.05;
            [~, power_fft, frq_fft, STA]  = STA_analysis(P(i).raster, P(i).antennal_movement, STA_window, fs,start_stim, stop_stim);
            P(i).power_fft = power_fft;
            P(i).frq_fft = frq_fft;
            P(i).STA = STA;
%             P(i).ISI = isi_hist(P(i).raster(:,start_stim:stop_stim),fs);
%             figure;
%             plot(frq_fft, power_fft);
%             xlabel('Frequency (in hertz)');
%             ylabel('Response magnitude');
%             title(join(split(P(i).stim_name,'_')," "));
%             xlim([0 350]);
            
%             figure; plot(STA);
%             figur
           
        end
    end
    
    
    % end
    % Covariance analysis
    
    for i = 1:no_of_protocols
        if P(i).stim_type == "blwgn"
            window = 0.05;
            [ev1, ev2] = cov_analysis(P(i).raster, P(i).antennal_movement, window, fs, start_stim, stop_stim);
%             P(i).cov_matrix = cov_matrix;
            P(i).ev1 = ev1;
            P(i).ev2 = ev2;
        end
    end
%%    % Latency to spike in square or step wave stimulus
    
    for i = 1:no_of_protocols
        
        if (P(i).stim_type == "sqr") || (P(i).stim_type == "step")
            latencies = spike_latency(P(i).antennal_movement,P(i).raster,start_stim, stop_stim,fs, P(i).complete_trials);
            P(i).latencies = latencies;
        end
        
    end
    
    %     mean_spike_latency = mean(extractfield(P, 'latencies'));
    %     spike_latency_SD = std(extractfield(P, 'latencies'));
    %     spike_latency_var = var(extractfield(P, 'latencies'));
    
    
    
% end
%% Write structure to another structure

% meta_struct(struct_num).moth = P;
% % amp_sweep(struct_num) = P(1);
% % blwgn(struct_num) = P(2);
% dec_frq_chirp(struct_num) = P;
% inc_frq_chirp(struct_num) = P;
% % T_sqr_27_01(struct_num,:) = struct2table(P, 'AsArray', true);
% struct_num = struct_num+1;

for i = 1:length(P)
    try
        if P(i).stim_type == "frq" || P(i).stim_type == "dec"
            chirp_meta_struct(chirp_meta_num).moth = P(i);
            chirp_meta_num = chirp_meta_num+1;
            
        elseif P(i).stim_type == "amp"
            ampsweep_meta_struct(ampsweep_meta_num).moth = P(i);
            ampsweep_meta_num = ampsweep_meta_num+1;
            
        elseif P(i).stim_type == "blwgn" || P(i).stim_type == "blwgn2"
            blwgn_meta_struct(blwgn_meta_num).moth = P(i);
            blwgn_meta_num = blwgn_meta_num+1;
            
        elseif P(i).stim_type == "sin" || P(i).stim_type == "sum"
            sin_meta_struct(sin_meta_num).moth = P(i);
            sin_meta_num = sin_meta_num +1;
            
        elseif P(i).stim_type == "sqr" || P(i).stim_type == "step"
            sqr_meta_struct(sqr_meta_num).moth = P(i);
            sqr_meta_num = sqr_meta_num+1;
            
        elseif (P(i).stim_type == "impulse")
            impulse_meta_struct(impulse_meta_num).moth = P(i);
            impulse_meta_num = impulse_meta_num + 1;
            
        else
            continue;
            
        end
    catch
        disp("failed to add recording: " + string(LUT_intra_row_idx));
        continue;
    end
end

% end
