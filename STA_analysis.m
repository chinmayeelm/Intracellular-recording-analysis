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
    t_STA = linspace(-(window*1000),0,length(STA));%-100:0.1:0;
    figure(); plot(t_STA, STA); %hold on;
    title ('Spike triggered average');
    ylabel 'Antennal movement (mm)';
    xlabel 'time (ms)';
    [STA_freq, power_fft, frq_fft] = fft_stim(STA, fs, window*fs);
%     histogram(STA_freq);
        

end
