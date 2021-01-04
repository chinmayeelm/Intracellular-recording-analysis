function [latency, mean_latency] = spike_latency(stimulus,raster_data,ON_dur,OFF_dur,fs, trials)
%SPIKE_LATENCY Summary of this function goes here
%   Detailed explanation goes here
    stim = stimulus(:, OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    resp = raster_data(:, OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    stim_SD = std(stim(1,1:(OFF_dur)*fs));
    
%     [m,~] = size(stim);
    
    for i=1:trials
        
        stim_on = find(stim(i,:) >= stim_SD);
        
        stim_on_loc = stim_on(1);
        spike_locs = (find(resp(i,:)==1));
        if length(spike_locs)==0
            continue;
        end
        
        first_spike_loc = spike_locs(1);
        
        latency(i) = (first_spike_loc - stim_on_loc)/fs;
    end
        mean_latency = mean(latency);
end

