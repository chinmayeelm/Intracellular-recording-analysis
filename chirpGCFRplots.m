function chirpGCFRplots(P_chirp_sweep)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

P = P_chirp_sweep(end);
total_trial_dur = P(1).ON_dur+2*P(1).OFF_dur;
time = linspace(0,total_trial_dur,P(1).single_trial_length);
figure('Color', 'w');

for i=1:length(P)
    
    ax1 = subplot(3,1,1); plot(time, -P(i).mean_movement, 'k', 'LineWidth',1); hold on;
    ylabel('Antennal position (deg)', 'FontSize',11);
    % yyaxis right; plot(time, P(i).intendedStimulus(1,:), 'r', 'LineWidth', 0.5);
    % ylabel('Generated position stimulus (a.u)');
    % grid on;
    %     lgd = legend(["1 s", "2 s","0.5 s", "0.5 s"],"Location","northeast","NumColumns",1);
    %     title(lgd, "Ramp duration");
    
    title(replace([P(1).date P(1).filename], '_','-'));
    ax1.Box = 'off';
    ax1.XAxis.Visible = 'off';
    ax1.XAxis.FontSize = 12;
    ax1.YAxis.FontSize = 12;
    
    % ax2 = subplot(4,1,2); plot(time(2:end),diff(P(i).intendedStimulus(1,:)),'k', 'LineWidth', 1); hold on;
    % ylabel('Velocity (a.u)', 'FontSize',11);
    % ax2.Box = 'off';
    % ax2.XAxis.Visible = 'off';
    % grid on;
    % 
    % f = linspace(0,150,P.ON_dur*P.fs);
    % % ax3 = subplot(3,1,2); plot(time(5*P.fs+1:30*P.fs),f,'k', 'LineWidth', 1, 'LineStyle','--'); hold on;
    % ax3 = subplot(3,1,2); plot(time(P.OFF_dur*P.fs+1:P.ON_dur)*P.fs),f,'k', 'LineWidth', 1, 'LineStyle','--'); hold on;
    % ylabel('Frequency (Hz)', 'FontSize',14);
    % ax3.Box = 'off';
    % ax3.XAxis.FontSize = 12;
    % ax3.YAxis.FontSize = 12;
    % ax2.XAxis.Visible = 'off';
    
    ax4 = subplot(3,1,3); plot(time, P(i).avg_gcfr,'k', 'LineWidth',1); hold on;
    ylabel('Mean Firing rate (Hz)', 'FontSize',14);
    xlabel('Time (s)', 'FontSize',14);
    ax4.Box = 'off';
    % grid on;
    ax4.XAxis.FontSize = 12;
    ax4.YAxis.FontSize = 12;
    
    linkaxes([ax1 ax3 ax4], 'x');
    xlim([3 Inf]);
    
end
% legend(ax1, arrayfun(@(x) replace(x.stim_name, "_"," "), P), 'Location', 'best');
% legend(ax1, 'boxoff');

end

