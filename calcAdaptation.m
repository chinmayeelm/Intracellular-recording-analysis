function [adaptation_coeff, rsquare] = calcAdaptation(P, fit_name)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fs = P(1).fs;
tonicGCFR = [];
gcfrDur = 1;
for i = 1:length(P)
    stim_name = split(P(i).stim_name);
    delta_t = str2double(stim_name(2));
    
    if delta_t==0.2500
        gcfr = P(i).avg_gcfr;
        [~,startPt] = max(gcfr);
        % startPt = P(i).OFF_dur*fs + delta_t*fs + 1;
        stopPt = startPt + gcfrDur*fs;
        % stopPt = P(i).OFF_dur*P(i).fs + delta_t*P(i).fs;
        tonicGCFR(i,:) = P(i).avg_gcfr(startPt : stopPt);
    end
end

descentFR = tonicGCFR(1,:);
% descentFR = (descentFR - min(descentFR))/(max(descentFR) - min(descentFR));
t1=linspace(1e-4, length(descentFR)/fs, length(descentFR));
% plot(t1,descentFR); hold on;

% idx = randperm(numel(t1),500);

[xData, yData] = prepareCurveData( t1, descentFR );
% xData = xData(idx);
% yData = yData(idx);

% Set up fittype and options.
ft = fittype( fit_name);
% ft = fittype( 'exp1' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Robust = 'Bisquare';
% opts.StartPoint = [9397215.94371203 -3.68090399369325];

% Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
[fitresult, gof] = fit( xData, yData, ft );

adaptation_coeff = fitresult;
rsquare = gof.rsquare;

% Plot fit with data.
figure();
plot( fitresult, xData, yData ); hold on;
xlim([0 Inf]);
% ylim([0 1.5]);
% legend( h, 't1 vs. descentFR', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'time (s)');
ylabel( 'Norm. mean firing rate');

str = sprintf(" k = %0.3f \n rsquare = %0.3f", adaptation_coeff, rsquare);
% text(1,1.5, str);
% text(0.25,1.5, str);
text(0.35,0.8, str);
% text(0.4,0.7, str);
grid on

title(join([string(P(1).date) replace(P(1).filename, '_',' ') P(1).stim_name], ' ' ));
% filename = join([replace(string(P(1).date),".","-") P(1).filename 'adaptation'], '_');
% cd('F:\Work\Analysis outputs\ramp_adaptation');
% saveas(gcf, filename , 'png');

end