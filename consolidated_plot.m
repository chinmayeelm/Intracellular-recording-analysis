function fig_handle = consolidated_plot(time, filtered_data_bp, hes_data, stim_fb, fs)
    fig_handle = figure();
    
    [p,l] = findpeaks(100*filtered_data_bp, "MinPeakHeight",0.3*max(100*filtered_data_bp), "MinPeakDistance", 0.003*fs);
    [p1,l1] = findpeaks(100*filtered_data_bp, "MinPeakHeight",0.1*max(100*filtered_data_bp), "MinPeakDistance", 0.003*fs);
    p_small = setxor(p,p1);
    l_small = setxor(l,l1);
    A1 = subplot(2,1,1);plot(time, filtered_data_bp*100, 'k');
    ylim([-10 80]);
    hold on; plot(l/fs, p, '.', 'MarkerEdgeColor', 'r');
%     plot(l_small/fs, p_small, '.', 'MarkerEdgeColor', 'b');
    
    title('Response to antennal movement')
    ylabel('Voltage (mV)')


%     subplot(3,1,2); hold on; plot(time, hes_data, 'Color', [0.6, 0.2,0]);
%     ylabel('Antennal movement (mm)')
%     A2 = gca;

    A3 = subplot(2,1,2); hold on; plot(time, stim_fb,'Color', [0.2,0.3,0.49]);
%     subplot(3,1,3); hold on; plot(time, stim_fb,'k');
    ylabel('Indenter feedback (V)')
    xlabel('time (s)')


    linkaxes([A1 A3], 'x');

%     net_movement_antenna = max(hes_data) - min(hes_data)

end