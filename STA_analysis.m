function [STA_freq, power_fft, frq_fft,STA] = STA_analysis(raster_data, stimulus, window, fs)

    [m,~] = size(raster_data);
    STA_freq = [];
    all_spike_triggers = [];
    for i=1:m
        spike_locs = find(raster_data(i,:)==1);
        if isempty(spike_locs)
            continue;
        end
        for j=1:length(spike_locs)
            if (spike_locs(j)-window*fs)<= 0 
                continue;
            end
            spike_triggers(j,:) = stimulus(i,(spike_locs(j)-window*fs):spike_locs(j));
            spike_triggers(j,:) = spike_triggers(j,:)-mean(spike_triggers(j,:));
%             plot(STA); hold on;
            
        end
%         [r,lags] =  xcorr(spike_triggers);
%         figure;
%         stem(lags,r);
        all_spike_triggers = [all_spike_triggers; spike_triggers];

    end
    
%     [r,lags] =  xcorr(all_spike_triggers);
%     figure;
%     stem(lags,r);
    

    STA = mean(all_spike_triggers);
    t_STA = linspace(-(window*1000),0,length(STA));%-100:0.1:0;
    figure(); plot(t_STA, STA); %hold on;
    title ('Spike triggered average');
    ylabel 'Antennal movement (mm)';
    xlabel 'time (ms)';
    [power_fft, frq_fft] = fft_stim(STA, fs, window*fs);
%     histogram(STA_freq);
    
%     rng('default');
% 
%     [cidx, ctrs] = kmeans(all_spike_triggers(1:50,:),2,'dist','corr','rep',3,'disp','final');
%     figure
%     for c = 1:16
%         subplot(4,4,c);
%         plot(times,all_spike_triggers((cidx == c),:)');
%         axis tight
end




