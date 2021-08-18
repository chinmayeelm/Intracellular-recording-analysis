function fig_handle = consolidated_plot(time, filtered_data_bp, hes_data, stim_fb, fs)
    fig_handle = figure();
    
    [p,l] = findpeaks(100*filtered_data_bp, "MinPeakHeight",0.25*max(100*filtered_data_bp));
    [p1,l1] = findpeaks(100*filtered_data_bp, "MinPeakHeight",0.1*max(100*filtered_data_bp));
    p_small = setxor(p,p1);
    l_small = setxor(l,l1);
    subplot(3,1,1);plot(time, filtered_data_bp*100, 'k');
    hold on; plot(l/fs, p, '.', 'MarkerEdgeColor', 'r');
    plot(l_small/fs, p_small, '.', 'MarkerEdgeColor', 'b');
    
    A1 = gca;
    title('Response to antennal movement')
    ylabel('Voltage (mV)')


    subplot(3,1,2); hold on; plot(time, hes_data, 'Color', [0.6, 0.2,0]);
    ylabel('Antennal movement (mm)')
    A2 = gca;

    subplot(3,1,3); hold on; plot(time, stim_fb,'Color', [0.2,0.3,0.49]);
%     subplot(3,1,3); hold on; plot(time, stim_fb,'k');
    ylabel('Indenter feedback (V)')
    xlabel('time (s)')
    A3 = gca;

    linkaxes([A1 A2 A3], 'x');

%     net_movement_antenna = max(hes_data) - min(hes_data)

end