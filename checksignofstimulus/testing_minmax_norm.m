%% Testing min-max normalization

fs = 10000;
t1 = 0:(1/fs):1;
tau = 0.1;
alpha = t1.*exp(-t1./tau);
s = [alpha -alpha];
t = linspace(0,2,length(s));
s_norm = miscFuncs.minmaxNorm(s);

[peakVal, tpeak] = max(s);
[peakVal_norm, tpeak_norm] = max(s_norm);

figure(1);
subplot(2,1,1);
plot(t,s); hold on;
ylim([min(s) max(s)]);
yline(peakVal/exp(1), 'g-');
tau_ = find(s(tpeak+1:end) < peakVal/exp(1),1,"first" );
xline(t(tau_+tpeak), 'g-');
tau = tau_/fs
title('Original signal')
text(0.5,peakVal/2, ['tau =' num2str(tau)])


[xData, yData] = prepareCurveData( t(tpeak:end), s(tpeak:end));
[fitresult, gof] = fit( xData, yData, 'power1' );
adaptation_coeff = fitresult.b;
rsquare = gof.rsquare;
str = sprintf(" k = %0.3f \n rsquare = %0.3f", adaptation_coeff, rsquare);
figure();
plot( fitresult, xData, yData ); hold on;
text(0.5,peakVal/2, str);
title('Original signal')




figure(1);
subplot(2,1,2)
% yyaxis right; 
plot(t,s_norm, 'k');
yline(peakVal_norm/exp(1), 'r-');
title('min-max normalization');
tau_norm_ = find(s_norm(tpeak_norm+1:end) < peakVal_norm/exp(1),1,"first" );
xline(t(tau_norm_+tpeak_norm), 'r-');
tau_norm = tau_norm_/fs
text(0.5,peakVal_norm/2, ['tau =' num2str(tau_norm)])

[xData_norm, yData_norm] = prepareCurveData( t(tpeak_norm:end), s_norm(tpeak_norm:end));
[fitresult_norm, gof_norm] = fit( xData_norm, yData_norm, 'power1' );
adaptation_coeff_norm = fitresult_norm.b;
rsquare_norm = gof_norm.rsquare;
str_norm = sprintf(" k = %0.3f \n rsquare = %0.3f", adaptation_coeff_norm, rsquare_norm);
figure();
plot( fitresult_norm, xData_norm, yData_norm ); hold on;
text(0.5,peakVal_norm/2, str_norm);
title('min-max normalised')



