cmp = parula(10);
for i=1:10
    % diff_sd = max_sd - min_sd;
    figure(1);
%     hold on;
    %     scatter( meta_struct(i).max_stim_sd, meta_struct(i).jitter, 10,cmp(i,:), 'filled');
    scatterhistogram(meta_table.max_stim_sd{i}, meta_table.jitter{i}, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Max STD in stimulus patterns (m)');
    ylabel('Jitter (s)');
    % title('58');
    
    
    figure(2);
    hold on;
    %     colormap turbo(11);
%     histogram(meta_struct(i).jitter,'BinWidth', 1000/meta_struct(i).fs,'FaceColor',cmp(i,:));
    histogram(meta_table.jitter{i},'BinWidth', 1000/meta_table.fs(i),'FaceColor',cmp(i,:));
    xlabel('Jitter (ms)');
    ylabel('Occurances');
    xlim ([0 2]);
% end    
    %%
    figure(3);
%     hold on;
    %     scatter(meta_struct(i).slope, meta_struct(i).jitter, 10,cmp(i,:), 'filled');
    scatterhistogram(meta_table.slope, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Slope of spike triggering stimulus pattern (mm)');
    ylabel('Jitter (ms)');
    
    figure(4);
%     hold on;
    %     scatter(meta_table.pos_amp, meta_struct(i).jitter, 10,cmp(i,:), 'filled');
    scatterhistogram(meta_table.pos_amp, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Amplitude of spike triggering stimulus pattern (mm/sec)');
    ylabel('Jitter (ms)');
    
    %     figure(5);
    %     hold on;
    % %     scatter(meta_struct(i).stim_pattern_r, meta_struct(i).jitter,10,cmp(i,:), 'filled');
    % scatterhistogram(meta_struct(i).stim_pattern_r, meta_struct(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    %     xlabel('Correlation with STA');
    %     ylabel('Jitter (ms)');
    
    figure(6);
%     hold on;
    %     scatter(meta_struct(i).time_to_prev_spike, meta_struct(i).jitter, 10,cmp(i,:), 'filled');
    scatterhistogram(meta_table.time_to_prev_spike, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Time to previous spike (ms)');
    ylabel('Jitter (ms)');
    
    %     colormap turbo(11);
    figure(7);
%     hold on;
    %     scatter(meta_struct(i).baseline_FR, mean(meta_struct(i).jitter), 40,cmp(i,:), 'filled');
    scatterhistogram(meta_table.baseline_FR, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Baseline firing rate (Hz)');
    ylabel('Jitter (ms)');
    ylim ([0 2]);
    
    %     colormap turbo(11);
    %     figure(8);
    %     hold on;
    % %     scatter(meta_struct(i).time_to_prev_spike_column, meta_struct(i).jitter(2:end), 10,cmp(i,:), 'filled');
    % scatterhistogram(meta_struct(i).time_to_prev_spike_column, meta_struct(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    %     xlabel('Time to previous spiking event (ms)');
    %     ylabel('Jitter (ms)');
    %     ylim ([0 2]);
    
    %     figure(9);
    %     hold on;
    % %     scatter(meta_struct(i).stim_ev_r, meta_struct(i).jitter,10,cmp(i,:), 'filled');
    %     scatterhistogram(meta_struct(i).stim_ev_r, meta_struct(i).jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    %     xlabel('Correlation with Eigen vector');
    %     ylabel('Jitter (ms)');
    
    %     figure(10);
    %     scatter(meta_struct(i).stim_pattern_freq, meta_struct(i).jitter, 10,cmp(i,:), 'filled');
    %     xlabel('Max frequency in spike triggering stimulus pattern (Hz)');
    %     ylabel('Jitter (ms)');
    %     hold on;
    
% end
