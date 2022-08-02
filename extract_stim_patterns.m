function stim_patterns = extract_stim_patterns(stimulus,raster, window, fs)

    [m,~] = size(raster);
    stim_patterns = [];
    
    for i=1:m
        spike_locs = find(raster(i,:)==1);
        if isempty(spike_locs)
            continue;
        end
        
        
        for j=1:length(spike_locs)-window*fs-1
            if (spike_locs(j)-window*fs)< 0 
                continue;
            end
            
            if (spike_locs(j)-spike_locs(j-1)) > window*fs && (spike_locs(j+1) - spike_locs(j)) > window*fs
                
                spike_triggers = stimulus(i,(spike_locs(j)-window*fs):(spike_locs(j)+window*fs));
                stim_patterns = [stim_patterns; spike_triggers];
            
            end
        end
    end
end