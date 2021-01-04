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
        all_spike_triggers = [all_spike_triggers; spike_triggers];

    end

    STA = mean(all_spike_triggers);
    t_STA = -100:0.1:0;
    figure(2); plot(t_STA, STA); hold on;
    [STA_freq, power_fft, frq_fft] = fft_stim(STA, fs, window*fs);
%     histogram(STA_freq);
        

end
