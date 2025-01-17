
function T = getJitterMeas(stimulus, raster, stim_window, spike_window, fs)

% stim_window = 0.01*fs; % 10 ms stimulus window
% spike_window = fs*2e-3; % 2 ms (on either side of spike = 4 ms) %earlier 1 ms on either side

[m,n] = size(raster);
min_no_trials_required = round(m/2);

locs_ref = find(raster(1,:)==1);


jitter_sd = [];
max_stim_sd = [];
min_sd = [];
max_sd = [];
freq = [];
slope = [];
pos_amp = [];
stim_pattern_r = [];
time_to_prev_spike = [];
max_coeff_var = [];
fidelity = [];
stim_ev_r = [];
pos_amp_sd = [];
STA = [];
STA_SD = [];
high_jitter_stim_patterns = [];
low_jitter_stim_patterns = [];
mid_jitter_stim_patterns = [];
mean_spiking_time = [];
slope_sd = [];

stim_peak_SD = [];


for j = 1:length(locs_ref)

    %     if locs_ref(j) > ref+spike_window
    ref = locs_ref(j);

    if ref>stim_window && (ref+spike_window) < length(stimulus)

        rows_with_multiple_ones = find(sum(raster(:,ref-spike_window : ref+spike_window), 2) > 1);
        if ~isempty(rows_with_multiple_ones)
            disp(rows_with_multiple_ones);
        end
        [rows, locs] = find(raster(:,ref-spike_window : ref+spike_window)==1);
        mult_spike_idx = ismember(rows, rows_with_multiple_ones);
        locs(mult_spike_idx) = [];
        locs_trials = ref+locs-spike_window-1; % 11 to get the actual location from raster


        if length(locs_trials)<min_no_trials_required || std(locs_trials)==0
            % disp(length(locs_trials));
            % disp(std(locs_trials));
            continue;
        end

        SD = (std(locs));
        jitter_sd = [jitter_sd; SD*1000/fs]; % * 1000 to convert to ms



        %             fidelity = [fidelity; 100*(length(spike_rows_trials)/P.complete_trials)];

        prev_spike = [];

        for k=1:length(locs)
            prev_spike(k) = find(raster(rows(k), 1:locs_trials(k)-1), 1,"last");
        end

        isi = (locs_trials - prev_spike');
        time_to_prev_spike = [time_to_prev_spike; median(isi)*1000/fs]; %in ms
        
        locs_trials_valid = (locs_trials-stim_window)>0;
        stim_patterns = stimulus(rows(locs_trials_valid), (locs_trials(locs_trials_valid)-stim_window):locs_trials(locs_trials_valid)-1);
        %             [stim_freq, ~, ~] = get_fft(mean(stim_patterns,1), fs, length(stim_patterns));
        %             freq = [freq; stim_freq];
        
        % for iter = 1:size(stim_patterns,1)
        %     [~, stim_locs]  = findpeaks(stim_patterns(iter,:), 'NPeaks', 2);
        %     % figure(j); findpeaks(stim_patterns(iter,:)); hold on;            
        %     if ~isempty(stim_locs)
        %         stim_peak_loc(iter) = stim_locs(end);
        % 
        %     else
        %         [~, stim_locs]  = findpeaks(-stim_patterns(iter,:), 'NPeaks', 2);
        %         % findpeaks(stim_patterns(iter,:), 'NPeaks', 2);
        %         stim_peak_loc(iter) = stim_locs(end);
        %         % stim_peak_loc(iter) = nan;
        %     end
        % 
        % end

        % stim_peak_SD = [stim_peak_SD std(stim_peak_loc)*1000/fs];

        if SD<=0.2*1e-3*fs
            low_jitter_stim_patterns = [low_jitter_stim_patterns; stim_patterns];

        elseif SD>0.2*1e-3*fs && SD<=0.35*1e-3*fs
            mid_jitter_stim_patterns = [mid_jitter_stim_patterns; stim_patterns];

        else
            high_jitter_stim_patterns = [high_jitter_stim_patterns; stim_patterns];

        end
        %
        %                         STA = mean([low_jitter_stim_patterns; mid_jitter_stim_patterns; high_jitter_stim_patterns], 1);
        %                         STA_SD = std([low_jitter_stim_patterns; mid_jitter_stim_patterns; high_jitter_stim_patterns], [],1);

        [max_pos, max_pos_ind] = max(stim_patterns,[], 2);
        [min_pos, min_pos_ind] = min(stim_patterns,[], 2);
        %             amp = stim_patterns(max_pos_ind) - stim_patterns(min_pos_ind);
        amp = abs(max_pos - min_pos);
        pattern_slope = amp*fs./abs(max_pos_ind - min_pos_ind);
        slope = [slope; median(pattern_slope)];
        pos_amp = [pos_amp; median(amp)];
        pos_amp_sd = [pos_amp_sd; std(amp)];
        slope_sd = [slope_sd; std(pattern_slope)];



        stim_patterns_sd = std(stim_patterns,0,1);
        stim_patterns_mean = mean(stim_patterns,1);
        stim_patterns_coeff_var = stim_patterns_sd./stim_patterns_mean;
        % min_sd = [ min_sd; min(stim_patterns_sd)];
        % max_sd = [max_sd; max(stim_patterns_sd)];
        max_stim_sd = [max_stim_sd; max(stim_patterns_sd)];
        max_coeff_var = [max_coeff_var; max(stim_patterns_coeff_var)];


    end

end


T= table({jitter_sd}, {max_stim_sd}, {max_coeff_var},{slope}, {pos_amp}, {pos_amp_sd}, {slope_sd}, {time_to_prev_spike}, {low_jitter_stim_patterns}, {mid_jitter_stim_patterns}, {high_jitter_stim_patterns});
T.Properties.VariableNames = {'jitter_sd', 'max_stim_sd', 'max_coeff_var','slope', 'pos_amp', 'pos_amp_sd', 'slope_sd', 'time_to_prev_spike', 'low_jitter_stim_patterns', 'mid_jitter_stim_patterns', 'high_jitter_stim_patterns'};

end
