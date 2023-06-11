cmp = parula(33);
for i=1:33
    
    figure(1);
    hold on;
    ax1 = scatter( S(i).max_stim_sd, S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.max_stim_sd{i}, meta_table.jitter{i}, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Max variation in stimulus patterns (mm)');
    ylabel('Jitter (s)');
    
    
    figure(2);
    hold on;
    
    ax2 = histogram(S(i).jitter,'BinWidth', 1000/S(i).fs,'FaceColor',cmp(i,:));
    xlabel('Jitter (ms)');
    ylabel('Occurances');
    xlim ([0 4]);
    
    
    figure(3);
    hold on;
    ax3 = scatter(S(i).slope/max(S(i).slope), S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.slope, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Slope of spike triggering stimulus pattern (mm/s)');
    ylabel('Jitter (ms)');
    
    figure(4);
    hold on;
    ax4 = scatter(S(i).pos_amp/max(S(i).pos_amp), S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.pos_amp, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Amplitude of spike triggering stimulus pattern (mm)');
    ylabel('Jitter (ms)');
    
%     figure(5);
%     hold on;
%     scatter(S(i).stim_pattern_r, S(i).jitter,10,cmp(i,:), 'filled');
%     % scatterhistogram(S(i).stim_pattern_r, S(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with STA');
%     ylabel('Jitter (ms)');
    
    figure(6);
    hold on;
    ax5 = scatter(S(i).time_to_prev_spike, S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.time_to_prev_spike, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Time to previous spike(ms)');
    ylabel('Jitter (ms)');
    
       
    
%     figure(8);
%     hold on;
%     scatter(S(i).time_to_prev_spike_column, S(i).jitter(2:end), 10,cmp(i,:), 'filled');
% %     scatterhistogram(S(i).time_to_prev_spike_column, S(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Time to previous spiking event (ms)');
%     ylabel('Jitter (ms)');
%     %     ylim ([0 2]);
    
%     figure(9);
%     hold on;
%     scatter(S(i).stim_ev_r, S(i).jitter,10,cmp(i,:), 'filled');
%     %     scatterhistogram(S(i).stim_ev_r, S(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with Eigen vector');
%     ylabel('Jitter (ms)');
    
    figure(10);
    ax6 = scatter(S(i).stim_pattern_freq, S(i).jitter, 10,cmp(i,:), 'filled');
    xlabel('Max frequency in spike triggering stimulus pattern (Hz)');
    ylabel('Jitter (ms)');
    hold on;
    xlim([0 350]);
    
end
