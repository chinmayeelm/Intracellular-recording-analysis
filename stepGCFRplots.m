function stepGCFRplots(P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

total_trial_dur = P(1).ON_dur+2*P(1).OFF_dur;
time = linspace(0,total_trial_dur,P(1).single_trial_length);
fs = P(1).fs;
newColors = [0.1765    0.1804    0.5137 1
             0.1765    0.1804    0.5137 0.5
             0.9137    0.3059    0.1059 0.5
             0.9137    0.3059    0.1059 1];
stim_name = string(extractfield(P, 'stim_name'));
step_pos = str2double(extractAfter(stim_name, "amp_ "));
[~,idx] = sort(step_pos);
P = P(idx);

[b,a] = butter(3,4/(fs/2), 'low');

figure('Color', 'w', 'WindowState','normal');

for i=1:length(P)

    S = abs(diff(P(i).intendedStimulus(1,1:P(i).single_trial_length)));
%     pos = max(abs(P(i).intendedStimulus(1,1:P(i).single_trial_length)))*100/2;
% 
%     x = [];
%     for j=1:length(S)-1
%         if (S(j)==0 && S(j+1)>0) || (S(j)>0 && S(j+1)==0)
%             x = [x j];
%         end
%     end
%     x = reshape(x,2,[]);
%     x(3,:) = x(2,:);
%     x(4,:) = x(1,:);
%     x = x/P(i).fs;
% 
%     y_stim = repmat([-1.2 -1.2 1.2 1.2]', 1,size(x,2));
%     y_gcfr = repmat([0 0 150 150]', 1,size(x,2));

    
    
    ax1 = subplot(3,1,1); plot(time, -P(i).mean_movement, 'LineWidth',1.5, 'Color',newColors(i,:)); hold on;
%     patch(x,y_stim,'k', 'FaceAlpha' , 0.05, 'EdgeColor', 'k');
    %     lgd = legend(["1 s", "2 s","0.5 s", "0.5 s"],"Location","northeast","NumColumns",1);
    %     title(lgd, "Ramp duration");
    ylabel('Position (deg)', 'FontSize',14, 'Rotation',0);
%     yyaxis right; ax5 = plot(time, P(i).intendedStimulus(1,:), '--', 'LineWidth', 0.5);
%     ylabel('Generated position stimulus (a.u)');
%     title(replace([P(1).date P(1).filename], '_','-'));
    ax1.Box = 'off';
    ax1.XAxis.Visible = 'off';
    ax1.YAxis.FontSize = 12;
    %grid on
    
    velocity = diff(-P(i).mean_movement)*fs;
    vel_filtered = filtfilt(b,a,velocity);
    ax2 = subplot(3,1,2); plot(time(2:end),vel_filtered,'Color', newColors(i,:), 'LineWidth', 1.5); hold on;
    ylabel('Velocity (deg/s)', 'FontSize',14, 'Rotation',0);
    ax2.Box = 'off';
    ax2.XAxis.Visible = 'off';
    ax2.FontSize = 12;
    % ax3 = subplot(4,1,3); plot(time(3:end),diff(P(i).intendedStimulus(1,:),2), 'LineWidth', 1); hold on;
    % ylabel('Accelaration (a.u)', 'FontSize',11);
    % ax3.Box = 'off';
    % ax3.XAxis.Visible = 'off';
    %grid on;
    
    ax4 = subplot(3,1,3); plot(time, P(i).avg_gcfr, 'LineWidth',1.5, 'Color',newColors(i,:)); hold on;
%     patch(x,y_gcfr, 'k', 'FaceAlpha' , 0.05, 'EdgeColor', 'k');
    ylabel('Mean Firing rate (Hz)', 'FontSize',14, 'Rotation',0);
    xlabel('Time (s)', 'FontSize',12);
    ax4.Box = 'off';
    ax4.FontSize = 12;
    %grid on;
    
    
    linkaxes([ax1 ax2 ax4], 'x');
    xlim([3 Inf]);    

end

%Steady state firing rate Vs Position

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
ylabel('Mean steady state firing rate (Hz)', 'FontSize',14, 'Rotation',0);
xlabel('Position (deg)', 'FontSize',14);


%
ylim([0 inf]);

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

end

