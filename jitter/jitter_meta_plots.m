cmp = parula(32);
fontSize = 12;

for i=1:32
    
    fig1 = figure(1);
    hold on;
    ax1 = scatter( S(i).amp_sd, S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.max_stim_sd{i}, meta_table.jitter{i}, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('STD of stimulus amplitude (mm)');
    ylabel('Jitter (s)');
    a1.YAxis.FontSize = 12;
    a1.XAxis.FontSize = 12;
    
    fig2 = figure(2);
    hold on;
    
    ax2 = histogram(S(i).jitter,'BinWidth', 1000/S(i).fs,'FaceColor',cmp(i,:));
    xlabel('Jitter (ms)');
    ylabel('Occurances');
    xlim ([0 4]);
    a2.YAxis.FontSize = 12;
    a2.XAxis.FontSize = 12;
    
    
    fig3 = figure(3);
    hold on;
    ax3 = scatter(S(i).slope/max(S(i).slope), S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.slope, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Slope of spike triggering stimulus pattern (mm/s)');
    ylabel('Jitter (ms)');
    a3.YAxis.FontSize = 12;
    a3.XAxis.FontSize = 12;
    
    fig4 = figure(4);
    hold on;
    ax4 = scatter(S(i).pos_amp/max(S(i).pos_amp), S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.pos_amp, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Amplitude of spike triggering stimulus pattern (mm)');
    ylabel('Jitter (ms)');
    a4.YAxis.FontSize = 12;
    a4.XAxis.FontSize = 12;
%     figure(5);
%     hold on;
%     scatter(S(i).stim_pattern_r, S(i).jitter,10,cmp(i,:), 'filled');
%     % scatterhistogram(S(i).stim_pattern_r, S(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with STA');
%     ylabel('Jitter (ms)');
    
    fig6 = figure(6);
    hold on;
    ax6 = scatter(S(i).time_to_prev_spike, S(i).jitter, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.time_to_prev_spike, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Time to previous spike(ms)');
    ylabel('Jitter (ms)');
    a6.YAxis.FontSize = 12;
    a6.XAxis.FontSize = 12;
       
    
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
    
    fig10 = figure(10);
    ax10 = scatter(S(i).stim_pattern_freq, S(i).jitter, 10,cmp(i,:), 'filled');
    xlabel('Max frequency in spike triggering stimulus pattern (Hz)');
    ylabel('Jitter (ms)');
    hold on;
    xlim([0 350]);
    a10.YAxis.FontSize = 12;
    a10.XAxis.FontSize = 12;
end
