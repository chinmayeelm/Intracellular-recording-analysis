function [vel_sorted,amplitude_sorted, max_FR_sorted] = rampGCFRplots(P,c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fs = P(1).fs;


stim_name = string(extractfield(P, 'stim_name'));
ramp_dur = str2double(extractAfter(stim_name, "ramp "));
[~,idx] = sort(ramp_dur);
P = P(idx);
n = numel(idx);
% n = 4;

[b,a] = butter(3,4/(fs/2), 'low');

%FR Vs Velocity


max_FR_mat = [];
velocity_mat = [];
delta_t_mat = [];
accel = [];
t_peakFR_mat= [];

for i=1:length(P)

    [onLoc, offLoc] = miscFuncs.findSSbounds(P(i).mean_movement, 0.9, 10, P(i).fs);
    ssBounds = [offLoc-1.5*P(i).fs  offLoc-0.5*P(i).fs];
    start_point = P(i).OFF_dur*fs+1;
    baselineBounds = [(P(i).OFF_dur-1.5)*fs+1 (P(i).OFF_dur-0.5)*fs];
    P(i).mean_vel_filtered = filtfilt(b,a,(diff(P(i).mean_movement)).*fs);
    stim_name = split(P(i).stim_name);
    delta_t = str2double(stim_name(2));
  
    baselineFR = mean(P(i).avg_gcfr(baselineBounds(1):baselineBounds(2)),"all");
    pos_ref = mean(P(i).mean_movement(baselineBounds(1):baselineBounds(2)));
    % pos_stim = mean(P(i).mean_movement(start_point + delta_t*P(i).fs + 1*P(i).fs:start_point + delta_t*P(i).fs + 3*P(i).fs));
    pos_stim = abs(mean(P(i).mean_movement(ssBounds(1):ssBounds(2))));
    amplitude = pos_ref - pos_stim;
    P(i).amplitude = amplitude;

    
    rampEndIdx = start_point+delta_t*fs;
    [xData, yData] = prepareCurveData((start_point:rampEndIdx)/fs, abs(P(i).mean_movement(start_point:rampEndIdx)));
    [ramp_fitresult, ramp_gof] = fit(xData, yData, 'poly1');
    velocity_value =  abs(ramp_fitresult.p1);
    ramp_fit_rsq = ramp_gof.rsquare;

    % figure; plot(ramp_fitresult, (start_point:rampEndIdx)/fs, P(i).mean_movement(start_point:rampEndIdx));

    velocity_mat = [velocity_mat velocity_value];

    acceleration = diff(P(i).mean_vel_filtered)*fs;
    P(i).mean_acc_filtered = filtfilt(b, a, acceleration);

    max_FR  = max(P(i).gcfr(:,start_point:onLoc), [], 2);
    max_FR_mat = [max_FR_mat max_FR];

    [mean_max_FR, t_peakFR] = max(P(i).avg_gcfr(:,start_point:onLoc));
    t_peakFR_mat = [t_peakFR_mat t_peakFR/fs];

    % figure;
    % protocolPlot(P(i)); hold on;
    % xline((start_point+t_peakFR)/fs, 'k--');
    % yline(mean_max_FR, 'k--');


end

velocity_mat =  repmat(velocity_mat, [P(i).complete_trials,1]);

% fig1 = plotPosVelAccGCFR(P, n, 0, 0, 0);
% title(replace(join([string(P(1).date) P(1).filename], " "), "_"," "));
% savefigures(P(1), "traces", fig1, "fig", 'D:\Work\Figures for presentation\uncategorized');
% savefigures(P, "traces", fig1, 1);

[~, idx] = sort(velocity_mat(1,:), 'ascend');
vel_sorted = velocity_mat(:,idx);
vel = vel_sorted(1,:);
max_FR_sorted = max_FR_mat(:,idx);
t_peakFR_sorted = t_peakFR_mat(:,idx);

% if mean(P(end).avg_gcfr(start_point:onLoc)) > 0
%         figure(26); scatter(vel_sorted(end), t_peakFR_sorted(end), 'filled'); hold on;
% end

meanMaxFR = mean(max_FR_sorted, 1);
err = std(max_FR_sorted, [],1);

amps = [P.amplitude];
amplitude_sorted = amps(idx);


[xData, yData] = prepareCurveData( vel_sorted, max_FR_sorted);
[fit_power, gof_power] = fit(xData, yData, 'power1');
k=fit_power.b;
pBounds = predint(fit_power, xData, 0.95, 'functional', 'on');
% err = std(max_FR_sorted, [],1);
yfit = fit_power(xData);

if gof_power.rsquare >=0.8
    % figure('WindowState', 'minimized');
    loglog(xData,yData, 'Color', [c 0.5],'Marker','.', 'MarkerSize',10, 'LineStyle','none'); hold on;
    loglog(xData, yfit, 'Color', c);
    if fit_power.b >=0
        fill([unique(xData); flip(unique(xData))], [unique(pBounds(:,1)); flip(unique(pBounds(:,2)))],...
            c, 'FaceAlpha',0.2, 'EdgeColor','none');
    else
        fill([unique(xData); flip(unique(xData))], [flip(unique(pBounds(:,1))); (unique(pBounds(:,2)))],...
            'k', 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    % text(1,100, string(gof_power.rsquare));
    % y_pred= feval(mdl, cell2mat(T_vel_fr.velocity(irow,:)));
    % plot(cell2mat(T_vel_fr.velocity(irow,:)), y_pred); hold on;
    % errorbar(cell2mat(T_vel_fr.velocity(irow)), cell2mat(T_vel_fr.fr(irow)), err, 'ro',"MarkerSize",3,...
    % "MarkerEdgeColor","none","MarkerFaceColor",'k');

ax=gca;
ax.XAxis.Scale = "log";
ax.YAxis.Scale = "log";
legend('', 'Box', 'off');
xlim([0.1 5])
xticks([0.1 0.5 1 2 3 4 5])
ylim([25 250])
yticks([25 50 100 150 200 250])
box off

ylabel('Peak firing rate (Hz)');
xlabel('Angular velocity ({\circ}/s)');
% savefigures(P(1), "stats", gcf, "png", 'D:\Work\Figures for presentation\uncategorized');
end

end