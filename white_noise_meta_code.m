% This code can be used to parse through all the valid white noise
% recordings

cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
filepath = readlines('blwgnValidData.txt');
dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);

cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");
% load("wn_chirp_files.mat");
% nfiles = height(wn_chirp_files);
% T_isi = table();
% T_jitter = table();

newcolors = [0 0 0;
    0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.4660 0.6740 0.1880];

% isi_neuron = nan(nfiles,3e4);
% max_isi_len = nan(nfiles,1);

corr_range = [];

for irow = 15%:nfiles
    irow
    % close all;

    P = getStructP(dataDirectory(irow), filename(irow), [nan nan], 1);

    % P = getStructP(wn_chirp_files.dataDirectory(irow), wn_chirp_files.chirp_filename(irow), [nan nan], 1);
   

    % plot_data(P(1), "stimulus", "raster");
    startPt = P(1).OFF_dur*P(1).fs+1;
    stopPt = (P(1).OFF_dur + P(1).ON_dur)*P(1).fs;
    
    
    % T_STA.stim(irow) = {P(1).antennal_movement(1,startPt:stopPt)};
    % T_STA.baseline_noise(irow) = {P(1).antennal_movement(1,1:startPt)};

    % Jitter
    % T = getJitterMeas(P.antennal_movement(:,startPt:stopPt), P.raster(:,startPt:stopPt), 0.01*P.fs, 0.002*P.fs, P.fs);
    % T_jitter = [T_jitter; T];
    

    % Chirp
    % [f,t] = instfreq(P_chirp(end).mean_movement(startPt_chirp:stopPt_chirp), P_chirp(end).fs, 'FrequencyLimits', [0 P_chirp(end).max_chirp_frq]);
    % time = linspace(0, P_chirp(end).ON_dur, P_chirp(end).ON_dur*P_chirp(end).fs); 
    % f_fit = fit(t,f, 'poly3');
    % % figure; plot(f_fit,t,f);
    % inst_chirp_freq = f_fit(time);
    % % figure; plot(inst_chirp_freq, P_chirp(end).avg_gcfr(startPt_chirp:stopPt_chirp));
    % 
    % wn_chirp_files.chirp_freq(irow) = {inst_chirp_freq};
    % wn_chirp_files.chirp_gcfr(irow) = {P_chirp(end).avg_gcfr(startPt_chirp:stopPt_chirp)};
    % pause;
    % chirpGCFRplots(P_chirp)
    


    
    % plot_data(P, "stimulus", "raster", "gcfr");
    % title(cellID(irow));
    % pause;
    
    % Adaptation
    % within_trial_adapt=[];
    % meanFR_across_trials = [];
    % meanFR_across_trials = sum(P(1).raster(:, startPt:stopPt), 2)./P.ON_dur;
    % k=1;
    % for j=startPt:P(1).fs:stopPt-P(1).fs+1 
    %     within_trial_adapt(:,k) = sum(P(1).raster(:,j:j+P(1).fs-1), 2);
    %     k=k+1;
    % end
    % 
    % Trow = table({meanFR_across_trials}, {within_trial_adapt});
    % T_adapt_wn = [T_adapt_wn; Trow];
    
    % baseline_activity = reshape(P.gcfr(:,1:startPt-1),1,[]);
    % stim_activity = reshape(P.gcfr(:, startPt:stopPt), 1, []);
    % [h,p] = kstest2(baseline_activity, stim_activity, 'Alpha',0.01);
    

    %ISI
    % raster = P(1).raster(:,startPt:stopPt);
    % [rows, cols] = find(raster);
    % isi = [];
    % 
    % 
    % for i=1:size(raster,1)
    %     idx = rows==i;
    %     splocs = cols(idx);
    %     isi_row = diff(splocs')*1000/P(1).fs;
    %     isi = [isi isi_row]; % to convert to ms
    % end

    % max_isi_len(irow) = length(isi);
    % isi_neuron(irow,1:length(isi)) = isi;

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
    % STA_window = 0.05;
    % load(join([cellID(irow) "_STE" ".mat"],""), 'STE');
    % STE = STE(:,end-STA_window*P(1).fs+1:end);
    % % load(join([cellID(irow) "_PSE" ".mat"],""), 'PSE');
    % fileID = join([cellID(irow) "_PSE_100" ".h5"],"");
    % PSE = h5read(fileID, '/PSE');
    % STA  = STA_analysis( P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), STA_window, P(1).fs, STE, PSE);
    % pause;
    % 
    % title(gca, replace([dataDirectory(irow), filename(irow)], "_", " "));
    % savefigures(P(1), "STA", gcf, "png", 'D:\Work\Figures for presentation\white-noise');
    % time = linspace(-100,0,length(STA));
    % Trow = table(dataDirectory(irow), filename(irow), {STA},{time},P.fs);
    % T_STA = [T_STA; Trow];

    %STE
    stim_window = 0.1;
    % STE = getIsolatedSTE(P(1).intendedStimulus(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs);
    % nSTE = size(STE,1);
    % fileID = join([cellID(irow) "_STE_intended_stimulus" ".mat"],"");
    % save(fileID, "STE",'-mat');
    % STE=[];

    %PSE
    nSTE = length(find(P(1).raster(:,startPt:stopPt)==1));
    % stim_window = 0.04;
    fileID = join([cellID(irow) "_PSE_100_hes" ".h5"],"");
    if exist(fileID, "file")
        continue;
    end
    h5create(fileID, '/PSE', [nSTE*100,stim_window*P(1).fs],'Datatype','single', ...
          'ChunkSize',[50 80],'Deflate',9);
    k = 1;
    for i = 1:100
        PSE = getPSE(P(1).antennal_movement(:,startPt:stopPt), stim_window, P.fs, nSTE);
        h5write(fileID, '/PSE', PSE, [k 1], size(PSE), [1 1]);
        k = k+nSTE;
    end

    

    % save(fileID, "PSE",'-mat');
    % PSE=[];

    % STA and EV
    % stim_window = 0.05;
    % [ev, STA] = cov_analysis(P(1).antennal_movement, P(1).raster, stim_window, P(1).fs,startPt, stopPt);
    % Trow = table(dataDirectory(irow), filename(irow), {STA}, {ev});
    % T_STA_EV = [T_STA_EV; Trow];

    

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

% T_isi = table(dataDirectory, filename, isi_neuron);
% T_isi.isi_neuron(:,max(max_isi_len)+1:end)=[];
% T_isi.Properties.VariableNames = {'date','filename','isi'};
% cellID = replace(join([T_isi.date T_isi.filename], " "),"_"," ");
% cellID = extractBetween(cellID, "2022-"," T" | " blwgn");
% T_isi.cellID = cellID;