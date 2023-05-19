function [adaptation_coeff, rsquare] = calcAdaptation(P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
fs = P(1).fs;
tonicGCFR = [];
for i = 1:length(P)
    gcfr = P(i).avg_gcfr;
    [~,startPt] = max(gcfr);
    stopPt = startPt + 4*fs;
    tonicGCFR(i,:) = P(i).avg_gcfr(startPt : stopPt);
end

descentFR = tonicGCFR(length(P),:);
descentFR = descentFR/max(descentFR);
t1=linspace(0.0001,4,4*fs+1);
plot(t1,descentFR); hold on;

idx = randperm(numel(t1),500);

[xData, yData] = prepareCurveData( t1, descentFR );
xData = xData(idx);
yData = yData(idx);

% Set up fittype and options.
ft = fittype( 'power1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Robust = 'Bisquare';
% opts.StartPoint = [9397215.94371203 -3.68090399369325];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts )

adaptation_coeff = fitresult.b;
rsquare = gof.rsquare;

% Plot fit with data.
% figure();
% h = plot( fitresult, xData, yData );
% legend( h, 't1 vs. descentFR', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 't1', 'Interpreter', 'none' );
% ylabel( 'descentFR', 'Interpreter', 'none' );
% text(2,1, join(string([adaptation_coeff rsquare])));
% grid on

title(join([string(P(1).date) replace(P(1).filename, '_',' ')], ' ' ));

end