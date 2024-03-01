% This code can be used to parse through all the valid white noise
% recordings

cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
filepath = readlines('blwgnValidData.txt');
dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);

T_STA_EV = table();

newcolors = [0 0 0;
    0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.4660 0.6740 0.1880];


for irow = 1:nfiles
    irow
    close all;

    P = getStructP(dataDirectory(irow), filename(irow), [nan nan], 1);
   

    plot_data(P(1), "stimulus", "raster");
    startPt = P(1).OFF_dur*P(1).fs+1;
    stopPt = (P(1).OFF_dur + P(1).ON_dur)*P(1).fs;
    
    %Jitter and fidelity
    % stim_window = 0.02*P(1).fs; % 20 ms stimulus window
    % spike_window = P(1).fs*2e-3; % 2 ms (on either side of spike = 4 ms) %earlier 1 ms on either side
    % min_no_trials_required = round(0.5*(P(1).complete_trials));
    % 
    % jitter = jitter_SD(P(1).raster(:,startPt:stopPt), stim_window, spike_window, min_no_trials_required, P(1).fs);
    % spikeProb = fidelity(P(1).raster, stim_window, spike_window);
    % Trow = table(dataDirectory(irow), filename(irow), {jitter}, {spikeProb});
    % T_jitter = [T_jitter; Trow];

    % pause;
    % figure;
    % histogram(jitter, 'BinWidth', 0.1);
    % xlim([0 4]); % 4 ms
    % pause;

    % STA
    STA_window = 0.1;
    STA  = STA_analysis(P(1).raster, P(1).antennal_movement, STA_window, P(1).fs, startPt, stopPt);
    title(gca, replace([dataDirectory(irow), filename(irow)], "_", " "));
    savefigures(P(1), "STA", gcf, "png", 'D:\Work\Figures for presentation\white-noise');
    Trow = table(dataDirectory(irow), filename(irow), STA);
    pause;

    % Corrcoef
    % for itrial = 1:P.complete_trials
    % 
    % stim_hes_unfilt = P.hes_data_unfilt(itrial,:) - P.hes_data_unfilt(itrial,1);
    % [x,y] = butter(10, 300/(P(1).fs/2), "low");
    % hes_filt = filtfilt(x,y, stim_hes_unfilt);
    % baseline_hes = hes_filt(1:startPt-1);
    % stim_hes_filt = hes_filt(startPt : stopPt);
    % 
    % stim_ifb = P.stim_ifb(itrial,startPt:stopPt) - P.stim_ifb(itrial,1);
    % % stim_intended = P.intendedStimulus(itrial,startPt:stopPt);
    % [h,p] = kstest2(baseline_hes, stim_hes_filt);    
    % c = corrcoef(stim_ifb, stim_hes_filt);
    % 
    % Trow = table(dataDirectory(irow),filename(irow),h, p, c(2));
    % T = [T; Trow];
    % end

    % PSD
    % figure;
    % plot(stim_ifb); hold on;
    % yyaxis right; plot(stim_hes_filt);
    % 
    % window = P(1).fs;
    % overlap = 0.5*window;
    % 
    % figure;
    % pwelch(stim_intended, [],[],[], P(1).fs); hold on;
    % pwelch(stim_ifb, [],[],[], P(1).fs); 
    % pwelch(stim_hes_unfilt, [],[],[], P(1).fs); 
    % pwelch(stim_hes_filt, [],[],[], P(1).fs); 
    % legend('generated stimulus', 'force feedback', 'HES unfiltered', 'HES filtered');
    % title(replace([dataDirectory(irow), filename(irow)], "_", " "));
    % xlim([0 0.4]);
    % colororder(newcolors);
    % savefigures(P(1), "PSD", gcf, "png", 'D:\Work\Figures for presentation\white-noise');


end