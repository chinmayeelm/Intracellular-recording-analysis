
load('D:\Work\Code\Intracellular-recording-analysis\tables\T_STA_withIsolatedSpikes.mat');
window = 0.1; % STA window

for i=4%1:32
    fs = T_STA.fs(i);
    sta = T_STA.STA{i};
    t = linspace(-window*fs,0,length(sta));
    [wt_sta, f_wt_sta] = cwt(sta, "amor", fs, VoicesPerOctave=48, ...
        FrequencyLimits=[1 500]);
    
    p_wt_sta = mean(abs(wt_sta),2);
    % figure; plot(f_wt_sta, p_wt_sta);
    p_wt = nan(50,length(f_wt_sta));
    for j=1:50
        sta_shift = circshift(sta,randi(1000));
        [wt,f] = cwt(sta_shift, "amor", fs, VoicesPerOctave=48, FrequencyLimits=[1 500]);
        p_wt(j,:) = mean(abs(wt),2);
        % figure(1); plot(t, sta_shift);
        % figure(2); plot(f,p_wt(j,:));
        % pause;
    end

    p_wt_mean = mean(p_wt,1);
    figure; hold on; plot(f, p_wt_sta)
    plot(f, p_wt_mean, 'k--');
    plot(f,p_wt_sta-p_wt_mean')
    legend('STA', 'Mean scrambled STA', 'Difference');
    box off
    legend box off
    xlabel('Frequency (Hz)')
    ylabel('Magnitude')
    title(i)
    pause;
end
