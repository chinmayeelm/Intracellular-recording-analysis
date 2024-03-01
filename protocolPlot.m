function protocolPlot(P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tickLabelSize = 10;
labelFontSize = 12;
neuronId = replace(string(join([P.date P.filename])), "_"," ");
c = lines(length(P));
for i=1:length(P)

    % figure;
    % stim_name = split(P(i).stim_name);
    % delta_t = str2double(stim_name(2));
    % [ss_onLoc, ss_offLoc] = miscFuncs.findSSbounds(P.mean_movement, 0.9, 10, P.fs);
    % ssBounds = [ss_onLoc+0.5*P(i).fs ss_offLoc-0.5*P(i).fs];
    % ssBounds = [ss_offLoc-1.5*P(i).fs  ss_offLoc-0.5*P(i).fs];
%{
    % stimulus plot
    A4 = subplot(2,1,1);
    hold on;
    % patch([P(i).OFF_dur,P(i).OFF_dur,P(i).OFF_dur+delta_t,P(i).OFF_dur+delta_t],[-1,1,1,-1], [0.96,0.9,0.9], 'FaceAlpha', 0, 'EdgeColor', 'k');
    % plot(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, 'DisplayName', neuronId, 'LineWidth',1);
    plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement, 'DisplayName', neuronId, 'LineWidth',1);
    % xline(P(i).OFF_dur+delta_t, 'k--', 'End of Ramp', 'LabelOrientation','aligned','LabelHorizontalAlignment','left');
    % label = {'SS ON', 'SS OFF'};
    % xline([ss_onLoc ss_offLoc]/P.fs, 'r--', label, 'LabelOrientation','aligned','LabelHorizontalAlignment','right');
    % xline(ssBounds/P.fs, 'r--', label, 'LabelOrientation','aligned','LabelHorizontalAlignment','right');
    % sdfill(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, std(P(i).antennal_movement,[],1), c(i,:)); %[0.6, 0.2,0]
    A4.Box = 'off';
    A4.XAxis.Visible = 'off';
    A4.FontSize = tickLabelSize;
    % A4.LineWidth =1;
    ylabel('Angular position (deg)','FontSize', labelFontSize);
    % ylim([-1 1]);
    A4.FontName = 'Calibri';

    % GCFR
    A3 = subplot(2,1,2);
    hold on;
    % patch([P(i).OFF_dur,P(i).OFF_dur,P(i).OFF_dur+delta_t,P(i).OFF_dur+delta_t],[-10,200,200,-10], [0.96,0.9,0.9],  'FaceAlpha', 0, 'EdgeColor', 'k');
    % plot(P(i).time(1:P(i).single_trial_length), P(i).avg_gcfr, 'LineWidth', 1);
    plot(P(i).time(1:P(i).single_trial_length), P(i).gcfr, 'LineWidth', 1);
    % xline(P(i).OFF_dur+delta_t, 'k--', 'End of Ramp', 'LabelOrientation','aligned','LabelHorizontalAlignment','left');
    % xline(ssBounds/P.fs, 'r--', label, 'LabelOrientation','aligned','LabelHorizontalAlignment','right');
    % plot(P(i).time(1:P(i).single_trial_length), P(i).gcfr);
    % sdfill(P(i).time(1:P(i).single_trial_length), P(i).avg_gcfr, std(P(i).gcfr,[],1), c(i,:)); %[0.4660 0.6740 0.1880]
    % yline(mean(P(i).avg_gcfr(1:4*P(i).fs)), 'k--');
    % ylim([0 150])
    A3.Box = 'off';
    A3.XAxis.Visible = 'on';
    A3.FontSize = tickLabelSize;
    % A3.LineWidth =1;
    A3.FontName = 'Calibri';
    ylabel('Mean firing rate (Hz)','FontSize', labelFontSize);
    xlabel('Time(s)','FontSize', labelFontSize);
    A3.Title.String = replace(join([string(P(1).date) P(1).filename P(i).stim_name], " "), "_"," ");
    if P(i).stim_name == "frq_chirp"
        frq = [zeros(1,P(i).OFF_dur*P(i).fs) linspace(0, P(i).max_chirp_frq, P(i).ON_dur*P(i).fs) zeros(1,P(i).OFF_dur*P(i).fs) 0];

        yyaxis right; plot(P(i).time(1:P(i).single_trial_length), frq);
    end

    linkaxes([A3,A4], 'x');
    A3.XLim = [3 20];

    % savefigures(P(i), "traces_all_trials", gcf, "fig", 'D:\Work\Figures for presentation\ramp_protocols_separate');
    % close all;
%}
    figure('WindowState','minimized');
    %Stimulus plot
    A4 = subplot(2,1,1);
    hold on;
    sdfill(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, std(P(i).antennal_movement,[],1), c(i,:)); %[0.6, 0.2,0]
    A4.Box = 'off';
    A4.XAxis.Visible = 'off';
    A4.FontSize = tickLabelSize;
    ylabel('Angular position (deg)','FontSize', labelFontSize);
    A4.FontName = 'Calibri';

    % GCFR
    A3 = subplot(2,1,2);
    hold on;
    sdfill(P(i).time(1:P(i).single_trial_length), P(i).avg_gcfr, std(P(i).gcfr,[],1), c(i,:)); %[0.4660 0.6740 0.1880]
    A3.Box = 'off';
    A3.XAxis.Visible = 'on';
    A3.FontSize = tickLabelSize;
    A3.FontName = 'Calibri';
    ylabel('Mean firing rate (Hz)','FontSize', labelFontSize);
    xlabel('Time(s)','FontSize', labelFontSize);
    A3.Title.String = replace(join([string(P(1).date) P(1).filename P(i).stim_name], " "), "_"," ");
    linkaxes([A3,A4], 'x');
    A3.XLim = [3 20];

    savefigures(P(i), "traces_sd", gcf, "png", 'D:\Work\Figures for presentation\step_protocols_separate');
    close all;
end

end