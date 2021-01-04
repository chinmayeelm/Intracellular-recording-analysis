function [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(stimulus, raster_data, fs, ON_dur, OFF_dur)

    stim = stimulus(OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    resp = raster_data(OFF_dur*fs:(OFF_dur+ON_dur)*fs);
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
    I_spike_freq = [];
    II_spike_freq = [];
    III_spike_freq = [];
    
    stim_freq = 1 ./ (diff(locs) / fs);
    stim_freq = sgolayfilt(stim_freq, 1,51);
    stim_period = 1./stim_freq;

    
    I = 1; II = 1; III=1;
    for k= 2:length(locs)
         stim_clips = stim(locs(k-1):locs(k));% figure(); subplot(2,1,1); plot(stim_clips) ;
         resp_clips = resp(locs(k-1):locs(k)); %subplot(2,1,2); plot(resp_clips);
    
         ind = find(resp_clips == 1);
         
         
         if length(ind)>=1
            I_spike_phase(I) = get_spike_phase(ind(1), fs, stim_period(k-1));  
            I_spike_freq(I) = stim_freq(k-1);
            I=I+1;
         end

         if length(ind) >=2
             II_spike_phase(II) = get_spike_phase(ind(2), fs, stim_period(k-1));
             II_spike_freq(II) = stim_freq(k-1);
             II = II+1;
         end

         if length(ind)>=3
             III_spike_phase(III) = get_spike_phase(ind(3), fs, stim_period(k-1));
             III_spike_freq(III) = stim_freq(k-1);
             III = III +1;
         end
    end
    
end

function spike_phase_val = get_spike_phase(index, fs, period)

    spike_phase_val = (index/fs)*((2*pi)/period);

end
    
    