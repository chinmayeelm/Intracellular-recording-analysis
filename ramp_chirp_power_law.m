% Comparing velocity and frequency ranges

load('D:\Work\Code\Intracellular-recording-analysis\tables\T_ramp_chirp.mat');

dataDirectory = extractBefore(ramp_chirp_fId.ramp_fileID, '_');
filename_ramp = extractAfter(ramp_chirp_fId.ramp_fileID, "_");
filename_chirp = extractAfter(ramp_chirp_fId.chirp_fileID, "_");

for irow = 1:height(ramp_chirp_fId)
    
    % P_ramp = getStructP(dataDirectory(irow), filename_ramp(irow), [nan nan], 1);
    P_chirp_all = getStructP(dataDirectory(irow), filename_chirp(irow), [nan nan], 1);

    P_chirp = P_chirp_all(end);
    
    % [velocity, amp, maxFR, k, rsq] = rampGCFRplots(P_ramp);
    
    fs = P_chirp.fs;

    startPt = P_chirp.OFF_dur*P_chirp.fs+1;
    stopPt = (P_chirp.OFF_dur + P_chirp.ON_dur)*P_chirp.fs;
    resp = P_chirp.raster(:,startPt:stopPt);
    stim = P_chirp.antennal_movement(1,startPt:stopPt);
    stim_ifb = P_chirp.stim_ifb(1,startPt:stopPt);
    ramp_chirp_fId.chirp_stim(irow) = {stim};
    % f_fit = miscFuncs.getFreqFit(stim_ifb, fs, P_chirp.ON_dur);
    % time = linspace(0,P_chirp.ON_dur,length(stim));
    % freq = f_fit(time); %P_chirp.inc_frq_chirp_f;
    % 
    % ramp_chirp_fId.velocity(irow) = {velocity};
    % ramp_chirp_fId.maxFR(irow) = {maxFR};
    % ramp_chirp_fId.k(irow) = k;
    % ramp_chirp_fId.rsq(irow) = rsq;
    % ramp_chirp_fId.freq(irow) = {freq};
    % ramp_chirp_fId.gcfr(irow) = {P_chirp.avg_gcfr(startPt:stopPt)};

end

%%

for irow = 1:height(ramp_chirp_fId)

    figure; 
    tiledlayout(3,1, 'TileSpacing','compact');

    ax1 = nexttile;
    plot(ramp_chirp_fId.freq{irow}, ramp_chirp_fId.chirp_stim{irow}, 'k');
    ylabel('Stimulus ({\circ})'); xlabel('Stimulus frequency (Hz)');
    axis padded;
    box off;
    set(gca().XAxis, 'Visible', 'off');

    ax2 = nexttile;
    plot(ramp_chirp_fId.freq{irow}, ramp_chirp_fId.gcfr{irow}, 'k'); 
    ylabel('Firing rate (Hz)'); xlabel('Stimulus frequency (Hz)');
    axis padded;
    box off;
    linkaxes([ax1 ax2], 'x');

    nexttile;
    [xData, yData] = prepareCurveData(cell2mat(ramp_chirp_fId.velocity(irow)), cell2mat(ramp_chirp_fId.maxFR(irow)));
    [f, gof] = fit(xData, yData, 'power1');
    yfit = f(xData);
    loglog(xData,yData, 'Color', [0 0 0 0.5],'Marker','.', 'MarkerSize',10, 'LineStyle','none'); hold on;
    loglog(xData, yfit, 'Color', 'k');
    text(mean(xData, 'all'), max(yData, [],'all') - 10, {gof.rsquare, f.b});
    ylabel('Peak firing rate (Hz)'); xlabel('Velocity ({\circ}/s)');
    axis padded;
    box off;

end

%% 

dataDirectory = extractBefore(ramp_ampsweep.ramp_fileID, '_');
filename_ramp = extractAfter(ramp_ampsweep.ramp_fileID, "_");
filename_chirp = extractAfter(ramp_ampsweep.chirp_fileID, "_");

for irow = 1:height(ramp_ampsweep)
    
    % P_ramp = getStructP(dataDirectory(irow), filename_ramp(irow), [nan nan], 1);
    P_ampsweep = getStructP(dataDirectory(irow), filename_chirp(irow), [nan nan], 1);

    P = P_ampsweep(1);
    
    % [velocity, amp, maxFR, k, rsq] = rampGCFRplots(P_ramp);
    
    fs = P.fs;

    startPt = P.OFF_dur*P.fs+1;
    stopPt = (P.OFF_dur + P.ON_dur)*P.fs;
    resp = P.avg_gcfr(:,startPt:stopPt);
    stim = P.mean_movement(1,startPt:stopPt);
    % stim_ifb = P.stim_ifb(1,startPt:stopPt);

    % figure;
    [p_stim, l_stim] = findpeaks(stim, 'MinPeakDistance', 0.1*P(1).fs);
    p_stim = 2*p_stim;
    window_size = 0.1*fs;
    p_resp = []; l_resp = [];

    for i = 1:length(p_stim)
        start_idx = max(1, l_stim(i) - window_size);
        end_idx = min(length(resp), l_stim(i) + window_size);
        [p, l] = max(resp(start_idx:end_idx));
        p_resp = [p_resp, p];
        l_resp = [l_resp, l + start_idx - 1];
    end
    % [p_resp, ~] = findpeaks(resp, 'MinPeakHeight', mean(resp)-10, 'MinPeakDistance', 0.1*P(1).fs, 'NPeaks',length(p_stim));
    [f,gof] = fit(p_stim', p_resp', 'power1');
    
    figure;
    t = linspace(0,P(1).ON_dur, length(stim));
    subplot(2,2,1); plot(t, stim, t(l_stim), p_stim/2, 'rx'); box off; axis padded;
    ylabel('Angular position (deg)');
    subplot(2,2,3); plot(t,resp, t(l_resp), p_resp, 'rx'); box off; axis padded;
    ylabel('Firing rate (Hz)'); xlabel('Time (s)');

    
    subplot(2,2,2);
    plot(f,p_stim, p_resp, 'o'); box off; axis padded;
    text(0.5, max(p_resp)-10, {gof.rsquare; f.b} );
    legend off;
    ylabel('Peak Firing rate (Hz)'); xlabel('Amplitude (deg)');
    
    subplot(2,2,4);
    [xData, yData] = prepareCurveData(cell2mat(ramp_ampsweep.velocity(irow)), cell2mat(ramp_ampsweep.maxFR(irow)));
    [f, gof] = fit(xData, yData, 'power1');
    yfit = f(xData);
    loglog(xData,yData, 'Color', [0 0 0 0.5],'Marker','.', 'MarkerSize',10, 'LineStyle','none'); hold on;
    loglog(xData, yfit, 'Color', 'k');
    text(mean(xData, 'all'), max(yData, [],'all') - 10, {gof.rsquare, f.b});
    ylabel('Peak firing rate (Hz)'); xlabel('Velocity ({\circ}/s)');
    axis padded;
    box off;

    % ramp_ampsweep.ampsweep_stim(irow) = {stim};
    % ramp_ampsweep.ampsweep_resp(irow) = {resp};
    ramp_ampsweep.amps(irow) = {p_stim};
    ramp_ampsweep.fr(irow) = {p_resp};
    ramp_ampsweep.ampsweep_k(irow) = f.b;
    ramp_ampsweep.ampsweep_rsq(irow) = gof.rsquare;



end

%% Chirp and amp sweep plots

% Chirp
figure; a1 = subplot(3,1,1); plot(P(end).time(1:P(end).single_trial_length), [P(end).antennal_movement]);
box off; axis padded;
ylabel('Antennal position (deg)');
[rr,ll] = find(P(end).raster == 1);
a2 = subplot(3,1,2); plot(ll/P(end).fs, rr, 'k|');
box off; axis padded;
ylabel('Trials');
a3 = subplot(3,1,3); sdfill(P(end).time(1:P(end).single_trial_length), P(end).avg_gcfr, std(P(end).gcfr, [], 1), 'k', "fill");
%plot(P(end).time(1:P(end).single_trial_length), P(end).avg_gcfr); 
box off; axis padded;
ylabel('Firing rate (Hz)');
xlabel('Time (s)');
linkaxes([a1,a2,a3], 'x');
% figure; plot(P(end).time(1:P(end).single_trial_length), P(end).mean_movement);
% hold on; plot(ll/P(end).fs, P(end).mean_movement(ll), 'r.', 'MarkerSize', 10);

figure;
for i=1:P(end).complete_trials
    ll = find(P(end).raster(i,:) == 1);
    positions = P(end).antennal_movement(i,ll);
    plot(P(end).time(1:P(end).single_trial_length), P(end).antennal_movement(i,:), 'LineWidth', 0.5);
    hold on;
    plot(ll/fs, positions, 'r.', 'MarkerSize', 20);
    box off; axis padded;
    ylabel('Antennal position (deg)');
    xlabel('Time (s)');
end

[positions,freq] = getFirstSpikeLocs(P(end).stim_ifb(:,startPt:stopPt), P(end).antennal_movement(:,startPt:stopPt),P(end).raster(:,startPt:stopPt), fs);
figure; plot(freq, positions, '.')
box off; axis padded;
xlabel('Frequency in the chirp signal (Hz)');
ylabel('Position (deg)');


%% amp sweep
figure; a1 = subplot(3,1,1); plot(P(1).time(1:P(1).single_trial_length), P(1).mean_movement);
box off; axis padded;
ylabel('Antennal position (deg)');

[rr,ll] = find(P(1).raster(1,:) == 1);
a2 = subplot(3,1,2); plot(ll/P(1).fs, rr, 'k|');
box off; axis padded;
ylabel('Trials');

a3 = subplot(3,1,3); sdfill(P(1).time(1:P(1).single_trial_length), P(1).avg_gcfr, std(P(1).gcfr, [], 1), 'k', "fill");
box off; axis padded;
ylabel('Firing rate (Hz)');
xlabel('Time (s)');
linkaxes([a1,a2,a3], 'x');
figure; plot(P(1).time(1:P(1).single_trial_length), P(1).antennal_movement(1,:));
hold on; plot(ll/P(1).fs, P(1).antennal_movement(1,ll), 'r.', 'MarkerSize', 10);

