function [ss_FR, position, fitresult, gof] = estimateTonicAdaptation(P, c)

[ss_onLoc, ss_offLoc] = miscFuncs.findSSbounds(abs(P.mean_movement), 0.9, 10, P.fs);

% ssBounds = [ss_onLoc+0.5*P.fs ss_onLoc+2*P.fs];
ssBounds = [ss_onLoc+0.1*P.fs ss_offLoc-0.5*P.fs];
gcfrTonic = P.avg_gcfr(ssBounds(1) : ssBounds(2));
tT = linspace(ssBounds(1)/P.fs, ssBounds(2)/P.fs, length(gcfrTonic));
t0 = linspace(0, length(gcfrTonic)/P.fs, length(gcfrTonic));

[xData, yData] = prepareCurveData( tT, gcfrTonic);

% Exponential fit
% mdl_exp = expFitPlots(xData, yData);
% y_pred_exp = feval(mdl_exp, t0);
% b = mdl_exp.Coefficients.Estimate(2);
% tau = 1/b;
% b_pval = mdl_exp.Coefficients.pValue(2);
% a_pval = mdl_exp.Coefficients.pValue(1);

[fitresult, gof] = expFitPlots(xData, yData);
y_pred_exp = feval(fitresult, tT);
tau = 1/fitresult.b;

c=[0,0,0]
t = linspace((1/P.fs), length(P.avg_gcfr)/P.fs, length(P.avg_gcfr));
% figure('WindowState','minimized');  
ax1= subplot(2,1,1);
plot(t, P.mean_movement, 'Color', c); hold on;
xline(ssBounds(1)/P.fs, 'k--');
ax1.XAxis.Visible = "off";
ylabel('Angular position ({\circ})');
box off;

ax2=subplot(2,1,2); 
plot(t, P.avg_gcfr, 'Color',c); hold on;
% ax=gca;
xline(ssBounds(1)/P.fs, 'k--');%,'Steady state begins','LabelHorizontalAlignment','right', 'LabelOrientation','horizontal');
plot(tT, gcfrTonic,	'g', tT, y_pred_exp, 'k--');

% plot(fitresult, tT, gcfrTonic, 'm');
% str1 = sprintf("{\\tau}_{1} (tonic adaptation) = %0.3f \np-value = %0.3f, %0.3f \nrsquare = %0.3f", abs(tau), a_pval, b_pval, mdl_exp.Rsquared.Ordinary);
str1 = sprintf("%0.3f", (tau));
text(tT(end)/2,gcfrTonic(round(length(gcfrTonic)/2)),str1);
ylabel("Mean firing rate (Hz)");
xlabel("Time (s)");
xlim([4 10]);
ylim([0 150])
box off;
linkaxes([ax1 ax2],'x');
% legend(['' '' str1], 'Box','off');
% savefigures(P, 'expfits', gcf, "png", 'D:\Work\Figures for presentation\expfits\step_toOffLoc');
% close all;


baselineBounds = [(P.OFF_dur-1.5)*P.fs+1 (P.OFF_dur-0.5)*P.fs];
pos_ref = mean(P.mean_movement(baselineBounds(1):baselineBounds(2)));
pos_stim = mean(P.mean_movement(ssBounds(1):ssBounds(2)));
position = pos_stim - pos_ref;
ss_FR = mean(P.avg_gcfr(ssBounds(1):ssBounds(2)));

end