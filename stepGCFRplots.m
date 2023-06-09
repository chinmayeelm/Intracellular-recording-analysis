function stepGCFRplots(P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

total_trial_dur = P(1).ON_dur+2*P(1).OFF_dur;
time = linspace(0,total_trial_dur,P(1).single_trial_length);

figure('Color', 'w', 'WindowState','maximized');
for i=1:length(P)
    
    ax1 = subplot(4,1,1); plot(time, -P(i).mean_movement, 'LineWidth',1); hold on;
    %     lgd = legend(["1 s", "2 s","0.5 s", "0.5 s"],"Location","northeast","NumColumns",1);
    %     title(lgd, "Ramp duration");
    ylabel('Antennal position (deg)', 'FontSize',11);
%     yyaxis right; ax5 = plot(time, P(i).intendedStimulus(1,:), '--', 'LineWidth', 0.5);
%     ylabel('Generated position stimulus (a.u)');
    title(replace([P(1).date P(1).filename], '_','-'));
    ax1.Box = 'off';
    ax1.XAxis.Visible = 'off';
    %grid on
    
    ax2 = subplot(4,1,2); plot(time(2:end),diff(P(i).intendedStimulus(1,:)), 'LineWidth', 1); hold on;
    ylabel('Velocity (a.u)', 'FontSize',11);
    ax2.Box = 'off';
    ax2.XAxis.Visible = 'off';
    %grid on;
    
    ax3 = subplot(4,1,3); plot(time(3:end),diff(P(i).intendedStimulus(1,:),2), 'LineWidth', 1); hold on;
    ylabel('Accelaration (a.u)', 'FontSize',11);
    ax3.Box = 'off';
    ax3.XAxis.Visible = 'off';
    %grid on;
    
    ax4 = subplot(4,1,4); plot(time, P(i).avg_gcfr, 'LineWidth',1); hold on;
    ylabel('Mean Firing rate (Hz)', 'FontSize',11);
    xlabel('Time (s)', 'FontSize',12);
    ax4.Box = 'off';
    %grid on;
    
    
    linkaxes([ax1 ax2 ax3 ax4], 'x');
    xlim([3 Inf]);    
end

%Steady state firing rate Vs Position
%{
figure;
steadystateFR = [];
position = [];
% FRsorted = [];
% position_sorted = [];
meanFR = [];
ssLoc = P(1).OFF_dur+3;

for i=1:length(P)
    stim_name = split(P(i).stim_name);
    delta_t = str2double(stim_name(2));
    pos_ref = mean(mean(P(i).antennal_movement(:,1:3*P(1).fs),2));
    pos_stim = mean(mean(P(i).antennal_movement(:, ssLoc*P(1).fs:(ssLoc+1)*P(1).fs),2));
    amplitude = pos_ref - pos_stim;
    position_values = repmat(amplitude, [P(i).complete_trials,1]);
    position = [position; position_values];
    
    steadystateFR = [steadystateFR; mean(P(i).gcfr(:,ssLoc*P(1).fs:(ssLoc+1)*P(1).fs),2)];
    
    
end

[position_sorted, idx] = sort(position, 'ascend');
FRsorted = steadystateFR(idx);
meanFR = median(reshape(FRsorted, P(i).complete_trials, []),1);
pos = unique(position_sorted);
% plot(pos, meanFR, '-o'); hold on;
% position_sorted = num2str(position_sorted, '%.3f');
boxchart(position_sorted, FRsorted, 'BoxWidth', 0.1, 'MarkerStyle','+'); hold on;
title(replace([P(1).date P(1).filename], '_','-'));
% boxplot(FRsorted, position_sorted);
ylabel('Firing rate (Hz)');
xlabel('Position (deg)');

FR =reshape(FRsorted, P(i).complete_trials, []);
if length(P) == 4
    plot(pos(1:2),median(meanFR(:,1:2),1), '-ro');
    plot(pos(3:4),median(meanFR(:,3:4),1), '-ro');
    scatter(position_sorted, FRsorted, 10,'k','filled','MarkerFaceAlpha',0.5);
    p1 = ranksum(FR(:,1), FR(:,2));
    p2 = ranksum(FR(:,3), FR(:,4));
else
    plot(pos(1:3),median(meanFR(:,1:3),1), '-ro');
    plot(pos(4:6),median(meanFR(:,4:6),1), '-ro');
    p1 = ranksum(FR(:,1), FR(:,3));
    p2 = ranksum(FR(:,4), FR(:,6));
end
text(-1,10,sprintf('%0.3f',p1));
text(0.5,10,sprintf('%0.3f',p2));
%
ylim([0 inf]);
%}

end

