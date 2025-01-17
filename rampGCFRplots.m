function [vel_sorted,amplitude_sorted, max_FR_sorted, trow] = rampGCFRplots(P,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    c = varargin{1};
else 
    c = [0 0 0];
end

fs = P(1).fs;


stim_name = string(extractfield(P, 'stim_name'));
ramp_dur = str2double(extractAfter(stim_name, "ramp "));
[~,idx] = sort(ramp_dur);
P = P(idx);
n = numel(idx);
% n = 4;

[b,a] = butter(3,4/(fs/2), 'low');

%FR Vs Velocity


max_FR_mat = nan(P(1).no_of_trials, length(P));
velocity_mat = nan(P(1).no_of_trials, length(P));
delta_t_mat = [];
accel = [];
t_peakFR_mat= [];

% figure(1);

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
    [xData, yData] = prepareCurveData((start_point:rampEndIdx)/fs, (P(i).mean_movement(start_point:rampEndIdx)));
    [ramp_fitresult, ramp_gof] = fit(xData, yData, 'poly1');
    velocity_value =  abs(ramp_fitresult.p1);
    ramp_fit_rsq = ramp_gof.rsquare;

    % subplot(2,1,1);
    % plot(ramp_fitresult, (start_point:rampEndIdx)/fs, P(i).mean_movement(start_point:rampEndIdx)); hold on;

    velocity_mat(1:P(i).complete_trials,i) = repmat(velocity_value,[P(i).complete_trials, 1]);

    acceleration = diff(P(i).mean_vel_filtered)*fs;
    P(i).mean_acc_filtered = filtfilt(b, a, acceleration);

    max_FR  = max(P(i).gcfr(:,start_point:onLoc), [], 2);
    max_FR_mat(1:length(max_FR),i) = max_FR;

    [mean_max_FR, t_peakFR] = max(P(i).avg_gcfr(:,start_point:onLoc));
    % t_peakFR_mat = [t_peakFR_mat t_peakFR/fs];

    % figure;
    % protocolPlot(P(i)); hold on;
    % xline((start_point+t_peakFR)/fs, 'k--');
    % yline(mean_max_FR, 'k--');


end



% figure();
% fig1 = plotPosVelAccGCFR(P, n, 0, 0, 0);
% title(replace(join([string(P(1).date) P(1).filename], " "), "_"," "));
% savefigures(P(1), "traces", gcf, "png", 'D:\Work\Figures for presentation\ramp');
% savefigures(P, "traces", fig1, 1);

[~, idx] = sort(velocity_mat(1,:), 'ascend');
vel_sorted = velocity_mat(:,idx);
vel = vel_sorted(1,:);
max_FR_sorted = max_FR_mat(:,idx);
% t_peakFR_sorted = t_peakFR_mat(:,idx);

% if mean(P(end).avg_gcfr(start_point:onLoc)) > 0
%         figure(26); scatter(vel_sorted(end), t_peakFR_sorted(end), 'filled'); hold on;
% end

meanMaxFR = mean(max_FR_sorted, 1);
err = std(max_FR_sorted, [],1);

amps = [P.amplitude];
amplitude_sorted = amps(idx);


[xData, yData] = prepareCurveData( vel_sorted, max_FR_sorted);

% semilog
% fittype_semilog = fittype('a + b*log10(x)',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a','b'});
% fitOpt = fitoptions('Method', 'NonlinearLeastSquares','Robust', 'Bisquare','StartPoint', [min(max_FR_sorted,[],"all") 10], 'Algorithm', 'Levenberg-Marquardt','MaxIter', 1000);
% [fit_semilog, gof_semilog] = fit(xData, yData, fittype_semilog, fitOpt);
% k = fit_semilog.b;
% rsq = gof_semilog.rsquare;
% trow = table(rsq, fit_semilog.b, fit_semilog.a);
% pBounds = predint(fit_semilog, xData, 0.95, 'functional', 'on');

% logistic

% fittype_logistic = fittype('a/(1+exp(-b.*(x-c)))',...
%     'dependent',{'y'},'independent',{'x'},...
%     'coefficients',{'a','b', 'c'});
% fitOpt = fitoptions('Method', 'NonlinearLeastSquares','Robust', 'Bisquare',...
%     'StartPoint', [mean(max_FR_sorted(:,end)) 1 0],...
%     'Algorithm', 'Levenberg-Marquardt','MaxIter', 10000);
% [fit_logistic, gof_logistic] = fit(xData, yData, fittype_logistic, fitOpt);
% rsq = gof_logistic.rsquare;
% k = fit_logistic.b;
% trow = table(rsq, fit_logistic.a, fit_logistic.b, fit_logistic.c);
% pBounds = predint(fit_logistic, xData, 0.95, 'functional', 'on'); hold on;

% Power

% fitOpt = fitoptions('Method', 'NonlinearLeastSquares','Robust', 'Bisquare','StartPoint', [1 0], 'Algorithm', 'Levenberg-Marquardt','MaxIter', 1000);
% [fit_power, gof_power] = fit(xData, yData, 'power1', fitOpt); %power1
% k=fit_power.b;
% rsq = gof_power.rsquare;
% trow = table(rsq, fit_power.a, fit_power.b);
% pBounds = predint(fit_power, xData, 0.95, 'functional', 'on'); hold on;

% linear

[fit_lin, gof_lin] = fit(xData, yData, 'poly1');
k = fit_lin.p1;
rsq = gof_lin.rsquare;
trow = table(rsq, fit_lin.p1, fit_lin.p2);
% pBounds = predint(fit_lin, xData, 0.95, 'functional', 'on'); hold on;


% err = std(max_FR_sorted, [],1);
% yfit = fit_lin(xData);


%{
neuronID = replace(join([P(1).date P(1).filename], " "), "_", " ");
% figure;
% if rsq >=0.8

    plot(xData,yData, 'Color', [c 0.5],'Marker','.', 'MarkerSize',10, 'LineStyle','none'); hold on;
    % plot(mean(vel_sorted,1),mean(max_FR_sorted,1), 'Color', [c 0.5],'Marker','.', 'MarkerSize',12, 'LineStyle','none','DisplayName', neuronID); hold on;
    f = plot(fit_logistic);
    set(f, 'Color', c);
    legend off;
    str = sprintf("Slope = %0.3f, rsquare = %0.3f", k, rsq);
    % plot(xData, yfit, 'Color', c, 'DisplayName', str); hold on;
    % fplot(fit_power, 'Color', c);
%     % if fit_power.b >=0
        % fill([unique(xData); flip(unique(xData))], [unique(pBounds(:,1)); flip(unique(pBounds(:,2)))],...
        %     c, 'FaceAlpha',0.2, 'EdgeColor','none');
%     % else
        fill([unique(xData); flip(unique(xData))], [(unique(pBounds(:,1))); flip(unique(pBounds(:,2)))],...
            'k', 'FaceAlpha',0.2, 'EdgeColor','none');
%     % end
    text(mean(xData),max(yData)-10, str);

% ax=gca;
% ax.XAxis.Scale = "log";
% ax.YAxis.Scale = "log";
legend('', 'Box', 'off');
% xlim([0.1 5])
% xticks([0.1 0.5 1 2 3 4 5])
ylim([0 250])
% yticks([25 50 100 150 200 250])
box off
axis padded;
ylabel('Peak firing rate (Hz)');
xlabel('Angular velocity ({\circ}/s)');
% end
%}
end