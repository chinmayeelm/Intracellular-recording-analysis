function [fitresult, gof] = expFitPlots(xData, yData)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


beta0_empirical_phasic = [max(yData), 0];
fo = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Robust', 'LAR', 'Algorithm', 'Trust-Region', 'MaxIter', 1000, 'StartPoint', beta0_empirical_phasic);
ft = fittype(@(a,b,x) a*exp(-b*x))
[fitresult, gof] = fit( xData, yData, ft, fo)
% figure; plot(fitresult, xData, yData, 'm');


% nlmfit
% model = @(b,x) (b(1)*exp(b(2)*x)+b(3));
% model = @(b,x) (b(1)*x.^b(2)); %power law
% beta0_empirical_phasic = [max(yData), 0, 0];
% options = optimset('MaxFunEvals', 10000);
% beta0_refined_phasic = fminsearch(@(beta) norm(yData - model(beta, xData)), beta0_empirical_phasic, options)
% statOpts = statset('MaxIter',10000);
% mdl_phasic = fitnlm(xData, yData, model, beta0_refined_phasic, "Options",statOpts);



end