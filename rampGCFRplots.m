function rampGCFRplots(P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
total_trial_dur = P(1).ON_dur+2*P(1).OFF_dur;
time = linspace(0,total_trial_dur,P(1).single_trial_length);
% colormap('winter');
stim_name = string(extractfield(P, 'stim_name'));
ramp_dur = str2double(extractAfter(stim_name, "ramp "));
[~,idx] = sort(ramp_dur);
P = P(idx);
n = numel(idx);

figure('Color', 'w', 'WindowState', 'normal');

switch n
    
    case 7
        newColors = [0.1765    0.1804    0.5137 1
                     0.1765    0.1804    0.5137 0.7
                     0.1765    0.1804    0.5137 0.5
                     0.1765    0.1804    0.5137 0.3
                     0.9137    0.3059    0.1059 0.5
                     0.9137    0.3059    0.1059 0.7
                     0.9137    0.3059    0.1059 1];
    case 6 
        newColors = [0.1765    0.1804    0.5137 1
                     0.1765    0.1804    0.5137 0.7
                     0.1765    0.1804    0.5137 0.5
                     0.9137    0.3059    0.1059 0.5
                     0.9137    0.3059    0.1059 0.7
                     0.9137    0.3059    0.1059 1];
    case 5
        newColors = [0.1765    0.1804    0.5137 1
                     0.1765    0.1804    0.5137 0.7
                     0.1765    0.1804    0.5137 0.5
                     0.9137    0.3059    0.1059 0.7
                     0.9137    0.3059    0.1059 1];
    otherwise 
        newColors = [0.1765    0.1804    0.5137 1
                     0.1765    0.1804    0.5137 0.5
                     0.9137    0.3059    0.1059 0.5
                     0.9137    0.3059    0.1059 1];
end
% colororder(newColors);
            
              
for i=1:length(P)
    
    ax1 = subplot(2,1,1); plot(time, -P(i).mean_movement, 'LineWidth',1, 'Color', newColors(i,:)); hold on;
    %     lgd = legend(["1 s", "2 s","0.5 s", "0.5 s"],"Location","northeast","NumColumns",1);
    %     title(lgd, "Ramp duration");
    ylabel('Antennal position (deg)', 'FontSize',12);
%     yyaxis right; ax5 = plot(time, P(i).intendedStimulus(1,:), '--', 'LineWidth', 0.5);
%     ylabel('Generated position stimulus (a.u)');
    title(replace([P(1).date P(1).filename], '_','-'));
    ax1.Box = 'off';
    ax1.XAxis.Visible = 'off';
%     colormap(winter) 
    %grid on;
    
%     ax2 = subplot(4,1,2); plot(time(2:end),diff(P(i).intendedStimulus(1,:)), 'LineWidth', 1); hold on;
%     ylabel('Velocity (a.u)', 'FontSize',11);
%     ax2.Box = 'off';
%     ax2.XAxis.Visible = 'off';
%     %grid on;
%     
%     ax3 = subplot(4,1,3); plot(time(3:end),diff(P(i).intendedStimulus(1,:),2), 'LineWidth', 1); hold on;
%     ylabel('Accelaration (a.u)', 'FontSize',11);
%     ax3.Box = 'off';
%     ax3.XAxis.Visible = 'off';
%     %grid on;
    
    ax4 = subplot(2,1,2); plot(time, P(i).avg_gcfr, 'LineWidth',1, 'Color', newColors(i,:)); hold on;
    ylabel('Mean Firing rate (Hz)', 'FontSize',12);
    xlabel('Time (s)', 'FontSize',12);
    ax4.Box = 'off';
%     colormap(winter) 
    
%     linkaxes([ax1 ax2 ax3 ax4], 'x');
    linkaxes([ax1 ax4], 'x');
    xlim([3 Inf]);
    
end


% legend(ax1, arrayfun(@(x) replace(x.stim_name, "_"," "), P), 'Location', 'best');
% legend(ax1, 'boxoff');
% filename = join([replace(string(P(1).date), ".", "-") P(1).filename], '_');
% cd('F:\Work\Analysis outputs\ramp_adaptation');
% saveas(gcf, filename, 'png');
%FR Vs Velocity
% figure;
max_FR = [];
velocity = [];
max_FR_sorted = [];
velocity_sorted = [];
% raster = P.raster(start_stim:stop_stim);
auc = zeros(1,length(P));
nSpikes = zeros(1,length(P));
avg_spikes = [];

for i=1:length(P)
    start_point = P(i).OFF_dur*P(i).fs+1;
    stim_end_point = (P(i).OFF_dur+P(i).ON_dur)*P(i).fs;
    stim_name = split(P(i).stim_name);
    delta_t = str2double(stim_name(2));
    baselineFR = mean(P(i).avg_gcfr(1:3*P(i).fs),"all");
    pos_ref = mean(-P(i).mean_movement(1:3*P(i).fs));
    pos_stim = mean(-P(i).mean_movement(start_point + delta_t*P(i).fs + 1*P(i).fs:start_point + delta_t*P(i).fs + 3*P(i).fs));
    amplitude = abs(pos_ref - pos_stim);
    velocity_values = repmat((amplitude/delta_t), [P(i).complete_trials,1]);
    velocity = [velocity; velocity_values];
    spikes = sum(P(i).raster(:,start_point:start_point+delta_t*P(i).fs),2)./delta_t;
    
    max_FR = [max_FR; max(P(i).gcfr, [], 2)];
    avg_spikes = [avg_spikes; spikes];
    
end

[velocity_sorted, idx] = sort(velocity, 'ascend');
max_FR_sorted = max_FR(idx);
avg_spikes_sorted = avg_spikes(idx);
meanMaxFR = median(reshape(max_FR_sorted, P(i).complete_trials, []),1);
meanSpikes = mean(reshape(avg_spikes_sorted, P(i).complete_trials, []),1);

vel = unique(velocity_sorted);
[vel_sorted,idx] = sort(unique(velocity), 'ascend');
auc = auc(idx);
nSpikes = nSpikes(idx);

% velocity_sorted = num2str(velocity_sorted, '%.3f');
figure;
boxchart(velocity_sorted, max_FR_sorted, 'BoxWidth', 0.03, 'MarkerStyle',"+"); hold on;
plot(vel, meanMaxFR);  hold on;
scatter(velocity_sorted, max_FR_sorted, 5,'k','filled','MarkerFaceAlpha',0.3);
% boxplot(max_FR_sorted, velocity_sorted);
% ylim([0 100]);
ylabel('Peak firing rate (Hz)', 'FontSize',12);
xlabel('Angular velocity (deg/s)', 'FontSize',12);
title(replace([P(1).date P(1).filename], '_','-'));

FR =reshape(max_FR_sorted, P(i).complete_trials, []);
p = ranksum(FR(:,1), FR(:,end))

end