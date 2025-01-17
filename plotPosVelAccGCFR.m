function figHandle = plotPosVelAccGCFR(P,ntrials, velFlag, accFlag, sdFillFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

switch ntrials

    case 7
        % newColors = [0.1765    0.1804    0.5137 1
        %     0.1765    0.1804    0.5137 0.7
        %     0.1765    0.1804    0.5137 0.5
        %     0.1765    0.1804    0.5137 0.3
        %     0.9137    0.3059    0.1059 0.5
        %     0.9137    0.3059    0.1059 0.7
        %     0.9137    0.3059    0.1059 1];

        newColors = [0.1765    0.1804    0.5137 1
            0.1765    0.1804    0.5137 0.7
            0.1765    0.1804    0.5137 0.5
            0    0.63    0.6 1
            0    0.63    0.6 0.7
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

        newColors = [0 0 0 1
            0 0 0 0.5
            0.8 0 0 0.5
            0.8 0 0 1];

        alphaVal = [1 0.85 0.70 0.55];
end

newColors = lines(length(P));

fs = P(1).fs;
total_dur = P(1).single_trial_length/fs;
time = linspace(0,total_dur,P(1).single_trial_length);

figHandle = gcf; % figure("WindowState","normal");
nSubplots = 2+velFlag+accFlag;


for iSubplot = 1:nSubplots
    ax(iSubplot) = subplot(nSubplots,1,iSubplot);
end

for iP=1:ntrials%length(P)

    
    iSubplot = 1;
    ax(iSubplot) = subplot(nSubplots,1,iSubplot); plot(time, P(iP).mean_movement, 'LineWidth',1, 'Color', newColors(iP,:));
    hold on;
    if sdFillFlag == 1
        sdfill(time, P(iP).mean_movement, std(P(iP).antennal_movement, [],1), newColors(iP,:));
    end
    iSubplot = iSubplot + 1;

    if velFlag == 1
        ax(iSubplot) = subplot(nSubplots,1,iSubplot); plot(time(2:end),P(iP).mean_vel_filtered,'Color', newColors(iP,:), 'LineWidth', 1); hold on;
        iSubplot = iSubplot + 1;
    end
   
    if accFlag == 1
        ax(iSubplot) = subplot(nSubplots,1,iSubplot); plot(time(3:end),P(iP).mean_acc_filtered,'Color', newColors(iP,:), 'LineWidth', 1); hold on;
        iSubplot = iSubplot + 1;
    end

    ax(iSubplot) = subplot(nSubplots,1,iSubplot); plot(time, P(iP).avg_gcfr, 'LineWidth',1, 'Color', newColors(iP,:)); hold on;
    if sdFillFlag == 1
        sdfill(time, P(iP).avg_gcfr, std(P(iP).gcfr, [],1), newColors(iP,:));
    end

end

ax(1).YLabel.String = 'Angular Position (deg)';
ax(1).Box = 'off';
ax(1).XAxis.Visible = 'off';
% title(replace(join([string(P(1).date) P(1).filename], " "), "_"," "));
ax(1).Title.String = replace(join([string(P(1).date) P(1).filename], " "), "_"," ");

if velFlag == 1
ax(2).YLabel.String = 'Angular Velocity (deg/s)';
ax(2).Box = 'off';
ax(2).XAxis.Visible = 'off';
end

if accFlag == 1
ax(iSubplot-1).YLabel.String = 'Angular accelaration (deg/s^{2})';
ax(iSubplot-1).Box = 'off';
ax(iSubplot-1).XAxis.Visible = 'off';
end

ax(iSubplot).YLabel.String = 'Mean Firing rate (Hz)';
ax(iSubplot).XLabel.String = 'Time (s)';
ax(iSubplot).Box = 'off';

linkaxes(ax, 'x');
xlim([3 Inf]);

end