cmp = parula(10);
for i=1:10
    
%     figure(1);
%     hold on;
%     scatter( blwgn_jitter_meta(i).max_stim_sd, blwgn_jitter_meta(i).fidelity, 10,cmp(i,:), 'filled');
%     %     scatterhistogram(meta_table.max_stim_sd{i}, meta_table.fidelity{i}, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Max STD in stimulus patterns (m)');
%     ylabel('Jitter (s)');
    
    
    
    figure(2);
    hold on;
    
    histogram(blwgn_jitter_meta(i).fidelity,'BinWidth', 1000/blwgn_jitter_meta(i).fs,'FaceColor',cmp(i,:));
    %     histogram(meta_table.fidelity{i},'BinWidth', 1000/meta_table.fs(i),'FaceColor',cmp(i,:));
    xlabel('Jitter (ms)');
    ylabel('Occurances');
    xlim ([0 2]);
    
    
    figure(3);
    hold on;
    scatter(blwgn_jitter_meta(i).slope, blwgn_jitter_meta(i).fidelity, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.slope, meta_table.fidelity, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Slope of spike triggering stimulus pattern (mm/s)');
    ylabel('Jitter (ms)');
    
    figure(4);
    hold on;
    scatter(blwgn_jitter_meta(i).pos_amp, blwgn_jitter_meta(i).fidelity, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.pos_amp, meta_table.fidelity, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Amplitude of spike triggering stimulus pattern (mm)');
    ylabel('Jitter (ms)');
    
%     figure(5);
%     hold on;
%     scatter(blwgn_jitter_meta(i).stim_pattern_r, blwgn_jitter_meta(i).fidelity,10,cmp(i,:), 'filled');
%     % scatterhistogram(blwgn_jitter_meta(i).stim_pattern_r, blwgn_jitter_meta(i).fidelity, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with STA');
%     ylabel('Jitter (ms)');
    
    figure(6);
    hold on;
    scatter(blwgn_jitter_meta(i).time_to_prev_spike, blwgn_jitter_meta(i).fidelity, 10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_table.time_to_prev_spike, meta_table.fidelity, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Time to previous spike (ms)');
    ylabel('Jitter (ms)');
    
       
    
%     figure(8);
%     hold on;
%     scatter(blwgn_jitter_meta(i).time_to_prev_spike_column, blwgn_jitter_meta(i).fidelity(2:end), 10,cmp(i,:), 'filled');
% %     scatterhistogram(blwgn_jitter_meta(i).time_to_prev_spike_column, blwgn_jitter_meta(i).fidelity, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Time to previous spiking event (ms)');
%     ylabel('Jitter (ms)');
%     %     ylim ([0 2]);
    
%     figure(9);
%     hold on;
%     scatter(blwgn_jitter_meta(i).stim_ev_r, blwgn_jitter_meta(i).fidelity,10,cmp(i,:), 'filled');
%     %     scatterhistogram(blwgn_jitter_meta(i).stim_ev_r, blwgn_jitter_meta(i).fidelity, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with Eigen vector');
%     ylabel('Jitter (ms)');
    
%     figure(10);
%     scatter(blwgn_jitter_meta(i).stim_pattern_freq, blwgn_jitter_meta(i).fidelity, 10,cmp(i,:), 'filled');
%     xlabel('Max frequency in spike triggering stimulus pattern (Hz)');
%     ylabel('Jitter (ms)');
%     hold on;
%     xlim([0 350]);
    
end
