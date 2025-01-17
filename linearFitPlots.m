function mdl_phasic = linearFitPlots(xData, yData)
%LINEARFITPLOTS Fits poly1 fit to phasic and tonic responses.
%   Detailed explanation goes here

% start_point = P.OFF_dur*P.fs+1;
% stim_name_parts = split(P.stim_name); 
% delta_t = str2double(stim_name_parts(2));
% rampEndIdx = start_point+delta_t*P.fs;
% 
% 
% 
% gcfr = P.avg_gcfr;%(start_point:start_point+rampEndIdx+4*P.fs);
% t = linspace((1/P.fs), length(gcfr)/P.fs, length(gcfr));
% figure('WindowState','minimized'); plot(t, gcfr); hold on;
% % xline([ss_onLoc ss_offLoc]/P.fs, 'k');
% xline(rampEndIdx/P.fs, 'k--','End of ramp','LabelHorizontalAlignment','right', 'LabelOrientation','horizontal');
% 
% 
% 
% [phasic_onLoc, ~] = miscFuncs.findSSbounds(P.avg_gcfr(start_point:rampEndIdx), 0.9, 10, P.fs);
% phasic_onLoc = phasic_onLoc + start_point;
% if ((rampEndIdx - phasic_onLoc) < 0.5*P.fs & delta_t >= 1) | isempty(phasic_onLoc) 
%     phasic_onLoc = rampEndIdx - delta_t*P.fs/2; 
% end
% phasic_offLoc = rampEndIdx - 0.01*P.fs; %10 ms before EoR
% 
% gcfrPhasic = (P.avg_gcfr(phasic_onLoc : phasic_offLoc));
% tPhasic = linspace((phasic_onLoc)/P.fs, phasic_offLoc/P.fs, length(gcfrPhasic));
% tP = linspace(0, length(gcfrPhasic)/P.fs, length(gcfrPhasic));

% [xData, yData] = prepareCurveData( tP, gcfrPhasic);
% [fitresult_phasic, gof_phasic] = fit( xData1, yData1, 'poly1');
% slope_phasic = fitresult_phasic.p1;

model = @(b,x) (b(1)*x + b(2));
beta0_empirical_phasic = [max(yData), 0.1];
options = optimset('MaxFunEvals', 1000);
beta0_refined_phasic = fminsearch(@(beta) norm(yData - model(beta, xData)), beta0_empirical_phasic, options);

mdl_phasic = fitnlm(xData, yData, model, beta0_refined_phasic);
% y_pred_phasic = feval(mdl_phasic, tP); 
% slope_phasic = mdl_phasic.Coefficients.Estimate(1);
% 
% plot(tPhasic, gcfrPhasic, tPhasic, y_pred_phasic, 'm');
% 
% 
% % plot(fitresult_phasic, xData1, yData1, 'm');
% str1 = sprintf("slope (phasic adaptation) = %0.3f \nrsquare = %0.3f", slope_phasic, mdl_phasic.Rsquared.Ordinary);
% legend(['' '' str1], 'Box','off');

%{
% Fit for decay in tonic response
% stepBounds = [rampEndIdx+1 rampEndIdx+2*P.fs];
stepBounds = [ss_onLoc+0.5*P.fs ss_offLoc-0.5*P.fs];
% xline(stepBounds/P.fs, 'r--')
gcfr_tonic = P.avg_gcfr(stepBounds(1) : stepBounds(2));
% tTonic = linspace((1/P.fs), length(gcfr_tonic)/P.fs, length(gcfr_tonic));
tTonic = linspace(stepBounds(1)/P.fs, stepBounds(2)/P.fs, length(gcfr_tonic));
[xData2, yData2] = prepareCurveData( tTonic, gcfr_tonic );
[fitresult_tonic, gof_tonic] = fit( xData2, yData2, 'poly1');
slope_tonic = fitresult_tonic.p1;
plot(fitresult_tonic, xData2, yData2, 'g');
str2 = sprintf("slope (tonic adaptation) = %0.3f \nrsquare = %0.3f", slope_tonic, gof_tonic.rsquare);
legend(['' '' str1 '' str2], 'Box','off');
%}
end