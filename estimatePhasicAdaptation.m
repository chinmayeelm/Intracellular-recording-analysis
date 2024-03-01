function [max_FR, velocity,fitresult_phasic, gof_phasic, ax2] = estimatePhasicAdaptation(P, c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

start_point = P.OFF_dur*P.fs+1;
P.stim_name
stim_name_parts = split(P.stim_name);
delta_t = str2double(stim_name_parts(2));
rampEndIdx = start_point+delta_t*P.fs;

if contains(P.stim_name, "step")
    [ss_onLoc, ss_offLoc] = miscFuncs.findSSbounds(abs(P.mean_movement), 0.9, 10, P.fs);
    rampEndIdx = ss_onLoc;
    delta_t = (rampEndIdx - start_point)/P.fs;
end

gcfr = P.avg_gcfr;
t = linspace((1/P.fs), length(gcfr)/P.fs, length(gcfr));


% Fit for decay in phasic response

[phasic_onLoc, ~] = miscFuncs.findSSbounds(P.avg_gcfr(start_point:rampEndIdx), 0.95, 10, P.fs);
phasic_onLoc = phasic_onLoc + start_point

if ((rampEndIdx - phasic_onLoc) < 0.5*P.fs & delta_t >= 1) | isempty(phasic_onLoc) 
    phasic_onLoc = rampEndIdx - delta_t*P.fs/2; 
end
phasic_offLoc = rampEndIdx - 0.01*P.fs; %10 ms before EoR

gcfrPhasic = (P.avg_gcfr(phasic_onLoc : phasic_offLoc));
tPhasic = linspace((phasic_onLoc)/P.fs, phasic_offLoc/P.fs, length(gcfrPhasic));
tP = linspace(1/P.fs, length(gcfrPhasic)/P.fs, length(gcfrPhasic));
max_FR = max(gcfrPhasic);
[xData, yData] = prepareCurveData( tP, gcfrPhasic);

% figure('WindowState','normal');

ax1 = subplot(2,1,1); plot(t, P.mean_movement, 'Color', c); hold on;
xline(rampEndIdx/P.fs, 'k--');
ax1.XAxis.Visible = "off";
ylabel('Angular position ({\circ})');
ylim([-1.2 0.2]);
box off

% Exponential fit
% mdl_exp = expFitPlots(xData, yData);
% y_pred_exp = feval(mdl_exp, tP); 
% b_phasic = mdl_exp.Coefficients.Estimate(2);
% tau_phasic = 1/b_phasic;
% b_pval_phasic = mdl_exp.Coefficients.pValue(2);
% a_pval_phasic = mdl_exp.Coefficients.pValue(1);

[fitresult_phasic, gof_phasic] = expFitPlots(xData, yData);
y_pred_exp = feval(fitresult_phasic, tP);
tau_phasic = 1/fitresult_phasic.b;


ax2 = subplot(2,1,2); plot(t, gcfr,'Color', c);  hold on; 
xline(rampEndIdx/P.fs, 'k--','EoR','LabelHorizontalAlignment','right', 'LabelOrientation','horizontal');
plot(tPhasic, gcfrPhasic, 'g', tPhasic, y_pred_exp, 'm');
% str1 = sprintf("{\\tau}_{1} (phasic adaptation) = %0.3f \np-value = %0.3f, %0.3f \nrsquare = %0.3f", abs(tau_phasic), a_pval_phasic, b_pval_phasic, mdl_exp.Rsquared.Ordinary);
str1 = sprintf("{\\tau}_{1} = %0.3f", abs(tau_phasic));
% legend(['' '' str1], 'Box','off');
text(tPhasic(end),max(gcfrPhasic),str1);
ylabel("Mean firing rate (Hz)");
xlabel("Time (s)");
box off;

linkaxes([ax1 ax2], 'x');
xlim([4 10]);
% savefigures(P, 'expfits', gcf, "png", 'D:\Work\Figures for presentation\expfits\');



[velocity, ~] = getSlope(P.mean_movement(start_point:rampEndIdx), P.fs, 0);



end