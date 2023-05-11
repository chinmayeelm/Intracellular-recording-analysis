function STA = STA_analysis(raster, stim, window, fs,start_stim, stop_stim)
    
    raster_data = raster(:,start_stim : stop_stim);
    stimulus = -stim(:,start_stim : stop_stim);
    
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
%             spike_triggers(j,:) = spike_triggers(j,:)- mean(spike_triggers(j,:));
%             plot(STA); hold on;
            
        end

        all_spike_triggers = [all_spike_triggers; spike_triggers];

    end
    
    pattern_length = size(all_spike_triggers,2);
    NstimPriors = size(all_spike_triggers,1);
    stimulus_prior = zeros([NstimPriors pattern_length]);
    r = randi([1 (length(stimulus_model)-pattern_length)],1,NstimPriors);
    for j=1:NstimPriors
        stimulus_prior(j,:) = stimulus_model(randperm(size(stimulus_model,1),1),r(j):r(j)+pattern_length-1);
    end
    
    avg_stim = mean(stimulus_prior,1);

    STA = mean(all_spike_triggers-avg_stim,1);
%     plot(STA);
%     t_STA = linspace(-(window*1000),0,length(STA));%-100:0.1:0;
%     figure(); plot(t_STA, STA); %hold on;
%     title ('Spike triggered average');
%     ylabel 'Antennal movement (mm)';
%     xlabel 'time (ms)';    

end




