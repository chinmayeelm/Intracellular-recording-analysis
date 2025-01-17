function corr_val = predictFR(stimulus,actual_gcfr, f_pspike, fs, stim_window, linearFilter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
predStepWin = 1; 
time = linspace(1,2,2*fs);

figure;
ax1 = subplot(2,1,1); plot(time, stimulus(1,:)); hold on;
ylabel('Stimulus (deg)');

pspikePred = nan(1,length(time));
timestamps = [];

for iwin = 1:predStepWin:(length(stimulus)-stim_window*fs-1)
    
    test_stim = (stimulus(1,iwin:iwin+stim_window*fs-1));
    % plot(time((iwin:iwin+stim_window*fs-1)), test_stim); hold on;
    similarity = dot(test_stim, linearFilter)/(norm(test_stim)*norm(linearFilter));
    timestamps = [timestamps iwin+stim_window*fs];
    pspikePred(iwin+stim_window*fs) = f_pspike(similarity);
end
box off;
set(gca().XAxis, 'Visible', 'off');

L = fs/5;
sigma = 0.005;

alpha = ((L-1)/(2*sigma*fs));
gauss_win = gausswin(L, alpha);


pspikePred_filt = (1/sum(gauss_win))*conv(pspikePred,gauss_win,'same');
% ax2 = subplot(2,1,2);plot(time, pspikePred); hold on;

ax2 = subplot(2,1,2);plot(time, pspikePred_filt,'k'); 

ylabel('Firing rate');
hold on;
plot(time,actual_gcfr(1,:), 'Color','#D95319','LineWidth',1);
title('Predicted firing rate with projection on linear filter')
box off;
legend('Predicted', 'Actual');

legend('Location','best', 'Box','off')
xlabel('Time (s)');

linkaxes([ax1 ax2], 'x');

corr = corrcoef(pspikePred, actual_gcfr, 'Rows','complete')
corr_val = corr(1,2);
end