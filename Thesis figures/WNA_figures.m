%% Raster plot
dataDirectory = "2022.09.16";
filename = "M1_N1_blwgn";

P = getStructP(dataDirectory, filename,[nan nan],1);

startPt = P(end).OFF_dur*P(end).fs+1;
stopPt = (P(end).OFF_dur + P(end).ON_dur)*P(end).fs;

plot_data(P(1), "stimulus", "raster");

%% Jitter histogram
c = parula(height(T_jitter_fidelity));
figure; colororder(c);
for i = 1%1:height(T_jitter_fidelity)
    histogram(cell2mat(T_jitter_fidelity.jitter(i)), 'BinWidth', 0.1,...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'Normalization','probability',...
        'DisplayName',replace(join([T_jitter_fidelity.date(i) T_jitter_fidelity.filename(i)]," "), "_", " "));
    hold on;
end
set(gca, 'YLimitMethod', 'padded');
xlim([-0.25 4]);
box off; ylabel('Probability'); xlabel('Jitter (ms)');
median_jitter = median(cell2mat(T_jitter_fidelity.jitter));
std_jitter = std(cell2mat(T_jitter_fidelity.jitter));

text_str = sprintf("Median \\pm STD = %0.2f \\pm %0.2f ms", median_jitter, std_jitter);
text(2,0.4, text_str);

%% Fidelity histogram
c = parula(height(T_jitter_fidelity));
figure; colororder(c);
for i = 1:height(T_jitter_fidelity)
    histogram(cell2mat(T_jitter_fidelity.fidelity(i)), 'BinWidth', 0.05,...
        'FaceAlpha', 0.5, 'EdgeAlpha', 0.5, 'Normalization','probability',...
        'DisplayName',replace(join([T_jitter_fidelity.date(i) T_jitter_fidelity.filename(i)]," "), "_", " "));
    hold on;
end
set(gca, 'YLimitMethod', 'padded');
xlim([-0.05 1]);
box off; ylabel('Probability'); xlabel('Fidelity');
median_fidelity = median(cell2mat(T_jitter_fidelity.fidelity));
std_fidelity = std(cell2mat(T_jitter_fidelity.fidelity));

text_str = sprintf("Median \\pm STD = %0.2f \\pm %0.2f", median_fidelity, std_fidelity);
text(0.4,0.4, text_str);

%% STA and EV

dataDirectory = "2022-08-17";
filename = "M1_N1_blwgn";
cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");

P = getStructP(dataDirectory, filename,[nan nan],1, 2e4/5, 0.005);
% plot_data(P, "stimulus", "gcfr");
stim_window = 0.04;
startPt = P(1).OFF_dur*P(1).fs+1;
stopPt = (P(1).OFF_dur + P(1).ON_dur - 2)*P(1).fs; %2s of test set
mean_FR = mean(P.avg_gcfr(startPt:stopPt-2*P(1).fs));

% STE_long = getIsolatedSTE(P.antennal_movement(:,startPt:stopPt), P.raster(:, startPt:stopPt), stim_window, P.fs);
% STE = getSTE(P.antennal_movement(:,startPt:stopPt), P.raster(:, startPt:stopPt), stim_window, P.fs);
% PSE_long = getPSE(P.antennal_movement(:,startPt:stopPt), stim_window, P.fs, size(STE_long,1)*100);
%%
% fileID_STE = join([cellID "_STE" ".mat"],"");
% load(fileID_STE);
% STE = STE(:,end-stim_window*P(1).fs+1:end);
% 
% fileID_PSE = join([cellID "_PSE_100_hes" ".h5"],"");
% PSE = h5read(fileID_PSE, '/PSE');
% PSE(:,1:end-stim_window*P.fs) = [];
% PSE = [];
% STA = STA_analysis(P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs, STE, PSE);
[eVal, ev, sig_evec_orth, sta, avg_stim, stc, pc, diff_cov] = cov_analysis(P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs, STE_long, PSE_long);
% [eVal, ev, ev_orth, STA, avg_stim, stc, pc, diff_cov] = cov_analysis(P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs);
beep;
% STA = STA_analysis(P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs, STE, PSE);
% figure(1);
% t = linspace(-100,0,2000);
% plot(t, miscFuncs.minmaxNorm_Minus1ToPlus1(STA), 'LineWidth',2, 'DisplayName','Intended'); hold on;

% [eVal, ev, STA, stc, pc, diff_cov] = cov_analysis(P(1).intendedStimulus(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs);

%% Plots

% Eigen values

min_null = -0.7351;
max_null =  0.7213;
figure;
scatter((1:length(eVal)),sort(eVal, 'descend'), 10, "filled");
hold on;
yline([min_null max_null], '--', {'min' 'max'});
xlabel('Eigen value index');
ylabel('Eigen value');
box off;
axis padded;
ylim([-1.5 1.5]);

% 
% Eigen vectors
t = linspace(-stim_window*1000,0,length(sta));
figure; plot(t,ev(:,1), t, ev(:,2)); box off;
hold on; plot(t, sta, 'k--');
legend('Eigen vector 1', 'Eigen vector 2', 'STA');
legend('boxoff');

% STC
% Downsample the matrices to save vector graphics
figure;
h = heatmap((stc), 'Colormap', parula);
h.YDisplayData = flipud(h.YDisplayData);
h.GridVisible = "off";
Labels = linspace(-stim_window*1e3,0,length(stc));
CustomLabels = string(Labels);  
CustomLabels(mod(Labels,10) ~= 0) = " ";
h.XDisplayLabels = CustomLabels;
h.YDisplayLabels = flip(CustomLabels);
title("Spike Triggered Covariance");

% PC
figure;
h = heatmap(pc, 'Colormap', parula);
h.YDisplayData = flipud(h.YDisplayData);
h.GridVisible = "off";
Labels = linspace(-stim_window*1e3,0,length(pc));
CustomLabels = string(Labels);
CustomLabels(mod(Labels,10) ~= 0) = " ";
h.XDisplayLabels = CustomLabels;
h.YDisplayLabels = flip(CustomLabels);
title("Prior covariance");

% diff cov
figure;
h = heatmap(diff_cov, 'Colormap', parula);
h.YDisplayData = flipud(h.YDisplayData);
h.GridVisible = "off";
Labels = linspace(-stim_window*1e3,0,length(diff_cov));
CustomLabels = string(Labels);
CustomLabels(mod(Labels,10) ~= 0) = " ";
h.XDisplayLabels = CustomLabels;
h.YDisplayLabels = flip(CustomLabels);
title("STC-PC");

%%
stim_window = 0.01;
PSE = PSE_long(:,end-stim_window*P(1).fs+1:end);
STE = STE_long(:,end-stim_window*P(1).fs+1:end);
STA = STA(end-stim_window*P(1).fs+1:end);
ev1 = ev(end-stim_window*P(1).fs+1:end,1);
%%
f_pspike_STA = constructNLD(STE, PSE, ev1', mean_FR);
predictFR(P(1).antennal_movement(1,stopPt+1:stopPt+2*P.fs),P(1).gcfr(1,stopPt+1:stopPt+2*P.fs), f_pspike_STA, P.fs, stim_window, ev1);

%% DPFD 


cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
filepath = readlines('blwgnValidData.txt');
dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);
c = parula(nfiles);
cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");

for irow = 14:nfiles
    
    irow 
    P = getStructP(dataDirectory(irow), filename(irow),[nan nan],1, 2e4/5, 0.005);
    % plot_data(P, "stimulus", "gcfr");
    stim_window = 0.04;
    startPt = P(1).OFF_dur*P(1).fs+1;
    stopPt = (P(1).OFF_dur + P(1).ON_dur - 2)*P(1).fs; %2s of test set
    mean_FR = mean(P.avg_gcfr(startPt:stopPt));

    fileID_STE = join([cellID(irow) "_STE" ".mat"],"");
    load(fileID_STE);
    STE = STE(:,end-stim_window*P(1).fs+1:end);
    %
    fileID_PSE = join([cellID(irow) "_PSE_100" ".h5"],"");
    if exist(fileID_PSE, 'file')
        PSE = h5read(fileID_PSE, '/PSE');
    else
        fileID_PSE_mat = join([cellID(irow) "_PSE" ".mat"],"");
        PSE_struct = load(fileID_PSE_mat);
        PSE = PSE_struct.PSE; 
    end
    PSE(:,1:end-stim_window*P.fs) = [];
    
    STA = STA_analysis(P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs, STE, PSE);
    % [eVal, ev, STA, stc, pc, diff_cov] = cov_analysis(P(1).antennal_movement(:,startPt:stopPt), P(1).raster(:,startPt:stopPt), stim_window, P(1).fs, STE, PSE);
    
    %{
    fig = figure; 
    tl = tiledlayout(size(ev,2),2);
    title(tl, cellID(irow));
    for iev=1:size(ev,2)
        nexttile;
        plot(ev(:,iev));
        box off
        nexttile;
        plot(diff(ev(:,iev)), 	'Color', "#D95319");
        box off
    end
    savefigures(P, "ev_derivatives_1", fig, ".png", "D:\Work\Figures for presentation\white-noise\EV_derivatives");

    
    %}
    % close all
    [f_pspike_STA, gof]  = constructNLD(STE, PSE, STA, mean_FR);
    title(cellID(irow))
    savefigures(P, "NLD", gcf, ".png", "D:\Work\Figures for presentation\white-noise\NLD");
    close all;
    T(irow,:) = table(cellID(irow), f_pspike_STA.p1, f_pspike_STA.p2, f_pspike_STA.p3, f_pspike_STA.p4, gof.rsquare);
    
    % f = plot(f_pspike_STA); hold on;
    % set(f, 'Color', c(irow,:));
    % legend off;
    % corr_val = predictFR(P(1).antennal_movement(1,stopPt+1:stopPt+2*P.fs),P(1).gcfr(1,stopPt+1:stopPt+2*P.fs), f_pspike_STA, P.fs, stim_window, STA)
    
end
T.Properties.VariableNames = {'cellID', 'p1','p2','p3','p4', 'rsq'};
%% STA from table
figure; hold on;
c = parula(height(T_STA));
colororder(c);
cellID = replace(join([T_STA.date T_STA.filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");
parfor i = 1:height(T_STA)
    plot(cell2mat(T_STA.time(i)), cell2mat(T_STA.STA(i)),...
        'LineWidth',1,...
        'DisplayName', cellID(i)); hold on;
end

box off;
axis padded;
set(gca, "LineWidth",1, "FontName", "Calibri");
ylim([-0.4 0.5])
xlabel('Time before spike (ms)'); ylabel('Antennal movement ({\circ})');

axes('Position', [0.2 0.6 0.4 0.3])
parfor i = 1:height(T_STA)
    plot(cell2mat(T_STA.time(i)), cell2mat(T_STA.STA(i)),...
        'LineWidth',1,...
        'DisplayName', cellID(i)); hold on;
end
ylim([-0.4 0.5])

%% STA plots with positon velocity acceleration classification
pos_vel_ind = find(T_STA.pos ==1 & T_STA.vel==1);
plot(cell2mat(T_STA.time(pos_vel_ind(1),:)),cell2mat(T_STA.STA(pos_vel_ind,:)), 'LineWidth',1);
box off;
title('Position and velocity sensitive')
nexttile;
plot(cell2mat(T_STA.time(2,:)),cell2mat(T_STA.STA(2,:)), 'LineWidth',1);
box off;
title('Position and acceleration sensitive')
nexttile;
plot(cell2mat(T_STA.time(11,:)),cell2mat(T_STA.STA(11,:)), 'LineWidth',1);
box off;
title('Acceleration sensitive')
set(gca().XAxis, 'Visible','off')
ax1 = gca;
ax2 = gca;
ax3 = gca;
ax4 = gca;
linkaxes([ax1 ax2],'x')
linkaxes([ax3 ax4],'x')
xlim([-30 0])




%% PSD of STAs

figure; hold on;
c = parula(height(T_STA));
colororder(c);

for i = 1:height(T_STA)
    % sta = miscFuncs.minmaxNorm_Minus1ToPlus1(cell2mat(T_STA.STA(i)));
    sta = (cell2mat(T_STA.STA(i)));
    % sta = sta(end-0.05*T_STA.fs(i)+1:end);
    % [pxx, f] = pwelch(sta, floor(length(sta)/2), floor(length(sta)/3), [], T_STA.fs(i));
    % pwelch(sta, [],[],[], T_STA.fs(i));
    %     plot(f, 10*log10(pxx/max(pxx)), 'LineWidth',1); hold on;
    %     ylim([-40 0]);
    %     xlim([0 300]);
    %     grid on;
    %     box off
    %     title('Power density estimate of STAs');
    %     ylabel('Normalized Power spectral density (dB/Hz)');
    % xlabel('Frequency (Hz)');
    % figure;plot(T_STA.time{i,1}, T_STA.STA{i,1});
    % figure;
    [stim_freq, power_fft, frq_fft] = fft_stim(cell2mat(T_STA.STA(i)), T_STA.fs(i));

    % xlim([0 0.3]);
    xlim([0 300]);
    T_STA.stim_freq(i) = stim_freq;
    T_STA.auc(i) = sum(sta(end-0.04*T_STA.fs(i)+1:end));
    T_STA.power_fft(i) = {power_fft};
    T_STA.frq_fft(i) = {frq_fft};
    % pause

end


%% PSD with Multitaper power spectral density and continuous wavelet transform

load('D:\Work\Code\Intracellular-recording-analysis\tables\T_STA_withIsolatedSpikes.mat');
% T_STA = T_STA_stim;
c = parula(32);
window = 0.1; % 40 ms
for i=1:32

    fs = T_STA.fs(i);
    sta = T_STA.STA{i};

    
    [p_mt, f_mt] = pmtm(sta, [], [], fs);
    [~, p_fft, f_fft] = fft_stim(sta, fs);
    [wt_sta, f_wt_sta] = cwt(sta, "amor", fs, VoicesPerOctave=48, FrequencyLimits=[1 500]);
    % figure(1);
    % cwt(sta, "amor", fs);
    % title(i);

    p_wt_sta = mean(abs(wt_sta), 2);
    fb = cwtfilterbank("SignalLength", length(sta), "Wavelet", "amor", "VoicesPerOctave", 48,"SamplingFrequency", fs, "FrequencyLimits", [1 500]);
    [p_wt_sta_ts, f_wt_sta_ts] = timeSpectrum(fb, sta);

   
    % figure(2);
    figure;
    
    subplot(2,1,1);
    plot(f_fft, p_fft/max(p_fft), 'DisplayName','FFT'); hold on;
    plot(f_mt, p_mt/max(p_mt), 'DisplayName', 'Multitaper');
    plot(f_wt_sta, p_wt_sta/max(p_wt_sta), 'Color', 'k', 'LineStyle','-', 'DisplayName', 'Wavelet');
    plot( f_wt_sta_ts, p_wt_sta_ts/max(p_wt_sta_ts),'DisplayName','tspec');
    legend('Box', 'off');
    box off;
    ylabel('Magnitude');
    xlim([0 300]);
    hold off;
    
  
    % figure(3);

    subplot(2,1,2);
    t= linspace(-100,0, length(sta));
    plot(t,sta, 'Color','k');
    ylabel('Antennal movement (deg)'); xlabel("Time before spike (ms)"); %hold on;
    box off;
    axis padded;
    xlim([-100 0]);
    title(i);
    
    pause;
end


%% PSD with continuous wavelet transform

load('D:\Work\Code\Intracellular-recording-analysis\tables\T_STA_withIsolatedSpikes.mat');
% T_STA = T_STA_stim;
c = parula(32);
window = 0.1; % 40 ms
ref = [60,120,180];
for i=1:32

    fs = T_STA.fs(i);
    sta = T_STA.STA{i};

    [wt_sta, f_wt_sta] = cwt(sta, "amor", fs, VoicesPerOctave=48, FrequencyLimits=[1 500]);
    % figure(1);
    figure;
    
    cwt(sta, "amor", fs);
    title(i);

    p_wt_sta = mean(abs(wt_sta), 2);
    [p,idx] = max(p_wt_sta);
    T_STA.pwt(i) = {p_wt_sta};
    T_STA.tuning_freq(i) = f_wt_sta(idx);
    T_STA.fwt(i) = {f_wt_sta};
    % figure(2);
    % figure;

    [~,ind] = min(abs(f_wt_sta(idx)-ref));
    % subplot(3,1,ind);
    % subplot(3,1,T_STA.freq_grp(i));
    figure;
    subplot(2,1,1);
    plot(f_wt_sta, p_wt_sta/max(p_wt_sta), 'Color', c(i,:), 'LineStyle','-', 'DisplayName', 'Wavelet'); hold on;
    % legend('Box', 'off');
    ylabel('Magnitude');
    
    set(gca().XAxis, 'Visible', 'off');
    % legend('location', 'best');
    % legend box off;

    
    plot(f_wt_sta(idx), p/p, 'Marker','o',"MarkerFaceColor",c(i,:), 'MarkerEdgeColor','none');
    % hold off;
    box off;
    axis padded;
    xlim([-10 300]);
    ylabel('Normalized Power/frequency (dB/Hz)');
    xlabel('Frequency (Hz)');
    % title(i);

    % sta_40ms = miscFuncs.minmaxNorm_Minus1ToPlus1(sta(end-0.04*fs+1:end));
    % t_40ms = linspace(-40, 0, length(sta_40ms));
    t = linspace(-100,0, length(sta));
    % figure(3);
    %
    % subplot(3,1,ind); hold on;
    % subplot(2,1,2);
    subplot(2,1,2);
    plot(t,sta, 'Color',c(i,:));
    ylabel('Antennal movement (deg)'); xlabel("Time before spike (ms)"); hold on;
    set(gca().XAxis, 'Visible', 'off');

    box off;
    axis padded;
    xlim([-50 0]);
    % title(i);
    %
    %
    % pause;
end

%%

[cfs,frq] = cwt(sta,"amor", fs);
tms = linspace(-100,0, length(sta));
figure
subplot(2,1,1)
plot(tms,sta);
box off;
ylim([-0.2 0.2]);
% title("Signal and Scalogram")
% xlabel("Time (s)")
ylabel("Antennal position (deg)");
subplot(2,1,2)
surface(tms,frq,abs(cfs))
axis tight
shading flat
xlabel("Time (s)")
ylabel("Frequency (Hz)")
set(gca,"yscale","log")
%% Swarmchart of tuning frequencies

figure; colororder(hot); swarmchart(ones(1,32), T_STA.tuning_freq, 20, 'filled', 'CData', T_STA.tuning_freq, 'XJitter', 'density', 'XJitterWidth', 1, 'MarkerFaceAlpha', 0.8);
colormap(gca, 'turbo')
xlim([0 2])
set(gca().XAxis, 'Visible', 'off');
camroll(90)
set(gca, 'YDir', 'reverse');

%% PSD vs response to frequency chirp
% Not matching tuning [8, 10, 14]
% Matching tuning [3, 13, 15]
% load('D:\Work\Code\Intracellular-recording-analysis\tables\wn_chirp_responses.mat');
figure;
% k=1;
tiledlayout(8,2, "TileSpacing","compact");
for i= 1:16
    % ax(k) = subplot(3,1,k);
    nexttile;
    plot(wn_chirp_files.chirp_freq{i}, wn_chirp_files.chirp_gcfr{i}, 'Color', c(idx(i),:), 'DisplayName', wn_chirp_files.cellID(i));
    

    yyaxis right; plot(wn_chirp_files.fwt{i}, wn_chirp_files.pwt{i}, 'Color', c(idx(i),:), 'LineStyle', '--', 'LineWidth', 2); 
    % hold on; plot(f_wt(idx), p/p, 'Marker','o',"MarkerFaceColor",c(i,:), 'MarkerEdgeColor','none');
    set(gca().YAxis(2), 'Color', [0.5 0.5 0.5]);
    set(gca().XAxis, 'Visible','off');
    box off;
    axis padded;
    xlim([0 300]);
    % k=k+1;
end
set(gca().XAxis, 'Visible','on');
% linkaxes(ax, 'x');
% xlim([20 120]);
xlabel('Frequency (Hz)');
ylabel('Average firing rate (Hz)');
yyaxis right; ylabel('Magnitude');
%% ISI plots

% load('D:\Work\Code\Intracellular-recording-analysis\tables\T_isi.mat');
n = height(T_isi);
c= parula(n);
med_isi = nan([n,1]);
sd_isi = nan([n,1]);
min_isi = nan([n,1]);
for i=1:n

    isi = T_isi.isi(i,:);
    med_isi(i) = median(isi, "omitnan");
    sd_isi(i) = std(isi, "omitnan");
    min_isi(i) = min(isi,[], "omitnan");
end

T_isi.median_ISI = med_isi;
T_isi.sd_ISI = sd_isi;
T_isi.min_ISI = min_isi;

figure;
for i=1:n
    histogram(T_isi.isi(i,:), 'BinWidth', 4, 'FaceColor', c(i,:), 'FaceAlpha', 0.3, 'Normalization', 'probability'); hold on;
end
box off;
axis padded;
ylabel('Probability');
xlabel('Inter-spike interval (ms)');

axes('Position', [0.4 0.4 0.5 0.5])
for i=1:n
    histogram(T_isi.isi(i,:), 'BinWidth', 4, 'FaceColor', c(i,:), 'FaceAlpha', 0.3, 'Normalization', 'probability'); hold on;
end
box off;
xlim([-3 50]);
ylim([-0.03 0.6]);

str = sprintf(string(["Least ISI = %0.2f" "\nMedian {\\pm} SD = %0.2f {\\pm} %0.2f" "\nn = %d"]), ...
    min(T_isi.min_ISI, [], "all", "omitnan"), median(T_isi.isi, "all", "omitnan"), std(T_isi.isi,[], "all", "omitnan"), height(T_isi), ...
    'interpreter', 'tex');
text(25,0.4, str);
figure; hold on;
[ncell,isiLen] = size(T_isi.isi);
pseudocell = ones(1,isiLen);
c = parula(ncell);
colororder(c);

for i = 1:ncell
    swarmchart(i*pseudocell, T_isi.isi(i,:), 5, 'XJitter', 'density');
end
% xticks(1:ncell);
% xticklabels(T_isi.cellID);


%% Response to frequency chirp and white noise

% load('D:\Work\Code\Intracellular-recording-analysis\tables\wn_chirp_responses.mat');

for irow = 1:height(wn_chirp_files)
    figure;
    plot(wn_chirp_files.chirp_freq{irow,1}, wn_chirp_files.chirp_gcfr{irow,1}); hold on;
    yyaxis right; plot(wn_chirp_files.fwt{irow,1}, wn_chirp_files.pwt{irow,1});
    box off;
    xlim([0 150]);
    pause;
end


%% Adaptation plots

load('D:\Work\Code\Intracellular-recording-analysis\tables\T_adapt_wn.mat');

c = parula(32);
figure; hold on;
for i = 1:32
    meanFR = mean(T_adapt_wn.meanFR_within_trial{i},1);
    plot((1:length(meanFR)), meanFR*100/(meanFR(1)), 'Color', c(i,:), 'Marker', '.', 'MarkerSize', 6, 'MarkerEdgeColor', c(i,:));
end
axis padded
box off
ylim([50 120])
xticks([1 5 10 15 20])
xlabel("Trial number");
ylabel("Percentage of mean firing rate");

%% Comparison of PSDs of stimulus signals

dataDirectory = "2022-08-17";
filename = "M1_N1_blwgn";
cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");
P = getStructP(dataDirectory, filename,[nan nan],1, 2e4/5, 0.005);
startPt = P.OFF_dur*P.fs + 1;
stopPt = (P.OFF_dur + P.ON_dur)*P.fs;
generated_stim = P.intendedStimulus(1,startPt:stopPt);
ifb_stim = P.stim_ifb(1,startPt:stopPt);
hes_stim = P.antennal_movement(1,startPt:stopPt);

% generated_stim_norm = miscFuncs.minmaxNorm_Minus1ToPlus1(generated_stim);
% ifb_stim_norm = miscFuncs.minmaxNorm_Minus1ToPlus1(ifb_stim);
% hes_stim_norm = miscFuncs.minmaxNorm_Minus1ToPlus1(hes_stim);
%%

stim_length = length(generated_stim);
win_length = stim_length/10;
overlap = 0.8*win_length;
[pxx_gen, f_gen] = pwelch(generated_stim,win_length,overlap,[], P.fs);
[pxx_ifb, f_ifb] = pwelch(ifb_stim,win_length,overlap,[], P.fs);
[pxx_hes, f_hes] = pwelch(hes_stim,win_length,overlap,[], P.fs);

figure;
hold on;
plot(f_hes, 10*log10(pxx_hes/max(pxx_hes)), 'LineWidth',2);
plot(f_ifb, 10*log10(pxx_ifb/max(pxx_ifb)), 'LineWidth',2);
plot(f_gen, 10*log10(pxx_gen/max(pxx_gen)), 'LineWidth',2);
ylabel('Normalized Power spectral density (dB/Hz)');
xlabel('Frequency (Hz)');
xlim([0 200]);
set(gca, 'FontName', 'Calibri');

%% Mean firing rate histogram

cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
filepath = readlines('blwgnValidData.txt');
dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);
c = parula(nfiles);
cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");

for irow = 1:nfiles
    
    irow 
    P = getStructP(dataDirectory(irow), filename(irow),[nan nan],1, 2e4/5, 0.005);
    % plot_data(P, "stimulus", "gcfr");
    if P(1).complete_trials~=15
        continue;
    else
    stim_window = 0.04;
    startPt = P(1).OFF_dur*P(1).fs+1;
    stopPt = (P(1).OFF_dur + P(1).ON_dur)*P(1).fs; %2s of test set
    % figure;
    % plot(P(1).mean_movement(6.5*P(1).fs:7*P(1).fs));
    % box off;
    % figure;
    % plot(P(1).mean_movement(11.5*P(1).fs:12*P(1).fs));
    % box off;
    % pause;
    plot_data(P(1), "stimulus", "raster");
    xlim([6.5 7]); 
    pause;

    xlim([11.5 12]);
    pause;
    end
end    


