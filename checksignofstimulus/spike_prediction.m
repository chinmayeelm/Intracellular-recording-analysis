    raster_data = blwgn(2).raster(1,:);
    actual_stim = blwgn(2).antennal_movement(1,:);
    sta = blwgn(1).STA(0.04*fs:end);

    t = linspace(0,15,length(raster_data));
    stim2 = blwgn(2).antennal_movement(1,:);
    
    stim_sta_conv = conv(stim2, sta, 'valid');
    [stim_sta_corr, lags] = xcorr(stim2, sta);
    
    
    figure;
    a1 = subplot(3,1,1); plot(t, stim2);
    ylabel 'Actual antennal movement';
    a2 = subplot(3,1,2); plot(stim_sta_conv);
    xlabel 'time (s)';
    ylabel 'Convolution of Stimulus and STA'
    
    a3 = subplot(3,1,3); stem(stim_sta_corr, lags);
    ylabel 'R';
    xlabel 'lags';
   
    
%     linkaxes([a1, a2], 'x');
    
    