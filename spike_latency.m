function latency = spike_latency(stimulus,raster_data,start_stim, stop_stim,fs, trials)
%SPIKE_LATENCY Summary of this function goes here
%   Detailed explanation goes here
    stim = stimulus(:, start_stim : stop_stim);
    resp = raster_data(:, start_stim : stop_stim);
    stim_SD = std(stim(1,1:start_stim));
    
%     [m,~] = size(stim);
    
    for i=1:trials
        
        stim_on = find(stim(i,:) >= stim_SD);
        
        stim_on_loc = stim_on(1);
        spike_locs = (find(resp(i,:)==1));
        if isempty(spike_locs)
            continue;
        end
        
        first_spike_loc = spike_locs(1);
        
        latency(i) = (first_spike_loc - stim_on_loc)/fs;
    end
        mean_latency = mean(latency);
end

