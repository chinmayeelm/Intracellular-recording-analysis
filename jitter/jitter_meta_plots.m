cmp = parula(32);
fontSize = 12;

figure;
for i=1:32
    
    i
    [wlow_jitter, whigh_jitter, med_jitter] = miscFuncs.getWhisker(T_jitter.jitter_sd{i,1});
    [wlow_amp_sd, whigh_amp_sd, med_amp_sd] = miscFuncs.getWhisker(T_jitter.pos_amp_sd{i,1});
    [wlow_slope, whigh_slope, med_slope] = miscFuncs.getWhisker(T_jitter.slope{i,1}/max(T_jitter.slope{i,1}));
    [wlow_amp, whigh_amp, med_amp] = miscFuncs.getWhisker(T_jitter.pos_amp{i,1}/max(T_jitter.pos_amp{i,1}));
    [wlow_tpsp, whigh_tpsp, med_tpsp] = miscFuncs.getWhisker(T_jitter.time_to_prev_spike{i,1});

    

    subplot(2,2,1);hold on;
    
    scatter( T_jitter.pos_amp_sd{i,1}, T_jitter.jitter_sd{i,1}, 10,cmp(i,:), 'filled', 'MarkerFaceAlpha', 0.5);hold on;
    % errorbar(med_amp_sd, med_jitter, wlow_jitter, whigh_jitter, wlow_amp_sd, whigh_amp_sd, 'Marker','o', 'MarkerFaceColor',cmp(i,:), 'Color', cmp(i,:), 'LineStyle','none');
    xlabel('STD of stimulus amplitude ({\circ})');
    ylabel('Jitter (s)');
    axis padded;


    % fig2 = figure(2);
    % figure;
    % hold on;
    % scatter( T_jitter.max_coeff_var{i,1}, T_jitter.jitter_sd{i,1}, 10,cmp(i,:), 'filled');
    % % plot( T_jitter.pos_amp_sd{i,1}, T_jitter.jitter_sd{i,1}, 'Color', cmp(i,:),'Marker','.','MarkerSize',20);
    % %     scatterhistogram(meta_table.max_stim_sd{i}, meta_table.jitter{i}, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    % xlabel('Coefficient of variation');
    % ylabel('Jitter (s)');
    % ax2 = gca;
    % % ax2.XAxis.Scale = "log";
    % % ax2.YAxis.Scale = "log";
    % ax2.YAxis.FontSize = 12;
    % ax2.XAxis.FontSize = 12;
    
    % fig2 = figure(2);
    % hold on;
    % 
    % histogram(cell2mat(T_jitter.jitter_sd(i)),'BinWidth', 1000/T_jitter.fs(i),'FaceColor',cmp(i,:));
    % xlabel('Jitter (ms)');
    % ylabel('Occurances');
    % xlim ([0 4]);
    % ax2= gca;
    % ax2.YAxis.FontSize = 12;
    % ax2.XAxis.FontSize = 12;
    
    
    % fig3 = figure(3);
    % figure;
    subplot(2,2,3);
    hold on;
    scatter(T_jitter.slope{i,1}/max(T_jitter.slope{i,1}), T_jitter.jitter_sd{i,1}, 10,cmp(i,:), 'filled', 'MarkerFaceAlpha', 0.5); hold on;
    % errorbar(med_slope, med_jitter, wlow_jitter, whigh_jitter, wlow_slope, whigh_slope, 'Marker','o', 'MarkerFaceColor',cmp(i,:),  'Color', cmp(i,:), 'LineStyle','none');
    % scatterhistogram(T_jitter.slope{i,1}/max(T_jitter.slope{i,1}), T_jitter.jitter_sd{i,1}, 'Color', cmp(i,:),'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    %     scatterhistogram(meta_table.slope, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Slope ({\circ}/s)');
    ylabel('Jitter (ms)');
    axis padded;
    
    % fig4 = figure(4);
    % figure;
    subplot(2,2,2);
    hold on;
    scatter(T_jitter.pos_amp{i,1}/max(T_jitter.pos_amp{i,1}), T_jitter.jitter_sd{i,1}, 10,cmp(i,:), 'filled', 'MarkerFaceAlpha', 0.5); hold on;
    % errorbar(med_amp, med_jitter, wlow_jitter, whigh_jitter, wlow_amp, whigh_amp, 'Marker','o', 'MarkerFaceColor',cmp(i,:),  'Color', cmp(i,:), 'LineStyle','none');
    % scatterhistogram(T_jitter.pos_amp{i,1}/max(T_jitter.pos_amp{i,1}), T_jitter.jitter_sd{i,1}, 'Color', cmp(i,:),'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    % ax4 = plot( T_jitter.pos_amp{i,1}/max(T_jitter.pos_amp{i,1}), T_jitter.jitter_sd{i,1}, 'Color', cmp(i,:),'Marker','.','MarkerSize',20);
    %     scatterhistogram(meta_table.pos_amp, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Amplitude ({\circ})');
    ylabel('Jitter (ms)');
    axis padded;
    
%     figure(5);
%     hold on;
%     scatter(T_jitter.stim_pattern_r, T_jitter.jitter,10,cmp(i,:), 'filled');
%     % scatterhistogram(T_jitter.stim_pattern_r, T_jitter.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with STA');
%     ylabel('Jitter (ms)');
    
    % fig6 = figure(6);
    % figure;
    subplot(2,2,4);
    hold on;
    scatter(T_jitter.time_to_prev_spike{i,1}, T_jitter.jitter_sd{i,1}, 10,cmp(i,:), 'filled', 'MarkerFaceAlpha', 0.5); hold on;
    % errorbar(med_tpsp, med_jitter, wlow_jitter, whigh_jitter, wlow_tpsp, whigh_tpsp, 'Marker','o', 'MarkerFaceColor',cmp(i,:),  'Color', cmp(i,:), 'LineStyle','none');
    % scatterhistogram(T_jitter.time_to_prev_spike{i,1}, T_jitter.jitter_sd{i,1}, 'Color', cmp(i,:),'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    % plot( T_jitter.time_to_prev_spike{i,1}, T_jitter.jitter_sd{i,1}, 'Color', cmp(i,:),'Marker','.','MarkerSize',20);
    %     scatterhistogram(meta_table.time_to_prev_spike, meta_table.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
    xlabel('Time to previous spike(ms)');
    ylabel('Jitter (ms)');
    axis padded;
       
    
%     figure(8);
%     hold on;
%     scatter(T_jitter.time_to_prev_spike_column, T_jitter.jitter(2:end), 10,cmp(i,:), 'filled');
% %     scatterhistogram(T_jitter.time_to_prev_spike_column, T_jitter.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Time to previous spiking event (ms)');
%     ylabel('Jitter (ms)');
%     %     ylim ([0 2]);
    
%     figure(9);
%     hold on;
%     scatter(T_jitter.stim_ev_r, T_jitter.jitter,10,cmp(i,:), 'filled');
%     %     scatterhistogram(T_jitter.stim_ev_r, T_jitter.jitter, 'HistogramDisplayStyle','smooth', 'LineStyle', '-', 'MarkerSize', 10, 'MarkerStyle', 'o', 'MarkerFilled', 'on' );
%     xlabel('Correlation with Eigen vector');
%     ylabel('Jitter (ms)');
    
    % fig10 = figure(10);
    % ax10 = scatter(T_jitter.stim_pattern_freq, T_jitter.jitter, 10,cmp(i,:), 'filled');
    % xlabel('Max frequency in spike triggering stimulus pattern (Hz)');
    % ylabel('Jitter (ms)');    
    % hold on;
    % xlim([0 350]);
    % a10.YAxis.FontSize = 12;
    % a10.XAxis.FontSize = 12;

    % pause;
end

%%
T_rho = table();

for i=1:32
    jitter_SD = T_jitter.jitter_sd{i,1};
    pos_amp_sd = T_jitter.pos_amp_sd{i,1};
    amp = T_jitter.pos_amp{i,1};
    slope = T_jitter.slope{i,1};
    time_to_prev_spike = T_jitter.time_to_prev_spike{i,1};

    T_rho.pos_amp_sd_jitter(i) =  corr(pos_amp_sd, jitter_SD, "type","Spearman");
    T_rho.amp_jitter(i) =  corr(amp, jitter_SD, "type","Spearman");
    T_rho.slope_jitter(i) =  corr(slope, jitter_SD, "type","Spearman");
    T_rho.time_to_prev_spike_jitter(i) =  corr(time_to_prev_spike, jitter_SD, "type","Spearman");

end

