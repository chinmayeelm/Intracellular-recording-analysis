function [I_spike_phase, II_spike_phase, III_spike_phase, I_spike, II_spike, III_spike] = spike_phase(stimulus, raster_data, fs, start_stim, stop_stim, stim_type)

stim = stimulus(start_stim : stop_stim);
resp = raster_data(start_stim : stop_stim);
mean_pos = mean(stim);

locs = [];
zc = [];
for j = 2:length(stim)
    if (stim(j-1)<=mean_pos && stim(j)>=mean_pos && stim(j+1)>mean_pos)
        zc(j) = 1;
    end
end
[val,locs ]= find(zc==1);

%     A1 =  subplot(2,1,1); plot(stim);hold on; plot(locs, stim(locs), 'rx'); yline(mean_pos); hold off;

stim_clips = [];
resp_clips = [];
I_spike_phase = [];
II_spike_phase = [];
III_spike_phase = [];
I_spike = [];
II_spike = [];
III_spike = [];




I = 1; II = 1; III=1;
for k= 2:length(locs)
    stim_clips = stim(locs(k-1):locs(k));
    resp_clips = resp(locs(k-1):locs(k));
    
    ind = find(resp_clips == 1);
    
    
    if stim_type == "frq"
        
        stim_freq = 1 ./ (diff(locs) / fs);
        stim_freq = sgolayfilt(stim_freq, 1,71);
        stim_period = 1./stim_freq;
        
        if length(ind)>=1
            I_spike_phase(I) = get_spike_phase(ind(1), fs, stim_period(k-1));
            I_spike(I) = stim_freq(k-1);
            I=I+1;
        end
        
        if length(ind) >=2
            II_spike_phase(II) = get_spike_phase(ind(2), fs, stim_period(k-1));
            II_spike(II) = stim_freq(k-1);
            II = II+1;
        end
        
        if length(ind)>=3
            III_spike_phase(III) = get_spike_phase(ind(3), fs, stim_period(k-1));
            III_spike(III) = stim_freq(k-1);
            III = III +1;
        end
        
    elseif stim_type == "amp"
        
        stim_freq = 1 ./ (diff(locs) / fs);
        stim_period = 1./stim_freq;
        
        if length(ind)>=1
            I_spike_phase(I) = get_spike_phase(ind(1), fs, stim_period(k-1));
            I_spike(I) = max(stim_clips) - min(stim_clips);
            I=I+1;
        end
        
        if length(ind) >=2
            II_spike_phase(II) = get_spike_phase(ind(2), fs, stim_period(k-1));
            II_spike(II) = max(stim_clips) - min(stim_clips);
            II = II+1;
        end
        
        if length(ind)>=3
            III_spike_phase(III) = get_spike_phase(ind(3), fs, stim_period(k-1));
            III_spike(III) = max(stim_clips) - min(stim_clips);
            III = III +1;
        end
    end
end

end

function spike_phase_val = get_spike_phase(index, fs, period)

spike_phase_val = mod((index/fs)*((2*pi)/period),2*pi);

end

