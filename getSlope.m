function [slope, rsq] = getSlope(signal,fs, plotFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    t = linspace(1/fs, length(signal)/fs, length(signal));
    
    [xData, yData] = prepareCurveData(t, signal);
    ft = 'poly1';
    [fitresult, gof] = fit(xData, yData, ft);
    slope = fitresult.p1;
    rsq = gof.rsquare;
    
    if plotFlag == 1
        figure;
        plot(fitresult, t, signal, 'k');
        str = sprintf("slope = %0.2f \n rsq = %0.2f", slope, rsq);
        legend(str, 'Location','best');
    end

end