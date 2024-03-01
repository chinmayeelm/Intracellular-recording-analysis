function [dist1,dist2] = distanceToFit(x, y, rampEndIdx, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[xData, yData] = prepareCurveData(x,y);
ft = 'exp2';
[fitresult, gof] = fit(xData, yData, ft)
figure; plot(fitresult, xData, yData); hold on;
tau1 = 1/fitresult.b;
tau2 = 1/fitresult.d;
rsq = gof.rsquare;

% str = sprintf("tau1 = %0.3f \n tau2 = %0.3f \n rsq = %0.3f", tau1, tau2, rsq);
% legend(str, 'Location','best');

f_tau1_curve_y = @(xData) fitresult.a*exp(fitresult.b*xData);
f_tau2_curve_y = @(xData) fitresult.c*exp(fitresult.d*xData);
tau1_curve_y = fitresult.a*exp(fitresult.b*xData);
tau2_curve_y = fitresult.c*exp(fitresult.d*xData);
fplot(f_tau1_curve_y, [-1 3], 'm', 'LineWidth',1);
fplot(f_tau2_curve_y, [-1 3], 'g', 'LineWidth',1);
start_value = 0.01;
intersection = rampEndIdx
% intersection = fzero(@(xData) f_tau1_curve_y(xData)-f_tau2_curve_y(xData), start_value)
xline(intersection, 'k--', 'intersection');
% tol = 3;
%
% if intersection <= 0.0 || intersection > x(end)
%     intersection1 = fzero(@(xData) tau1_curve_y(xData)-tol, start_value)
%     intersection2 = fzero(@(xData) tau2_curve_y(xData)-tol, start_value)
%     intersection = min([intersection1,intersection2]);
% intersection = rampEndIdx;
% end

dist1 = 0;
dist2 = 0;
dist1_for_curve2 = 0;
dist2_for_curve2 = 0;
i=1;
for i = 1:length(xData)

    if i<=intersection*fs
    dist1 = dist1 + (((yData(i) - tau1_curve_y(i))^2)/tau1_curve_y(i)^2);
    dist2 = dist2 + (((yData(i) - tau2_curve_y(i))^2)/tau2_curve_y(i)^2);
    % else
    %     dist1_for_curve2 = dist1_for_curve2 + (((yData(i) - tau1_curve_y(i))^2)/tau1_curve_y(i)^2);
    %     dist2_for_curve2 = dist2_for_curve2 + (((yData(i) - tau2_curve_y(i))^2)/tau2_curve_y(i)^2);
    % i=i+1;
    end
end


str = sprintf("intersection at %0.2f \ndistance to 1^{st} exponential =  %0.2f \ndistance to 2^{nd} exponential = %0.2f", intersection, dist1, dist2);
% str2 = sprintf("intersection at %0.2f \ndistance to 1^{st} exponential =  %0.2f \ndistance to 2^{nd} exponential = %0.2f", intersection, dist1_for_curve2, dist2_for_curve2);
text(1.5, max(yData)-30, str);

xline(rampEndIdx, 'r--', 'End of Ramp', 'LabelHorizontalAlignment','left');
legend({'data', 'Sum of 2 exponentials fit', 'First exponential', 'Second exponential', '', ''},'Box','off')
% title('Curve 2 distance')
% title(sprintf("Tol = %0.1f", tol));
end