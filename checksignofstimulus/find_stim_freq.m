function stim_frequencies = find_stim_freq(stimulus,ON_dur, OFF_dur, fs)
%FIND_STIM_FREQ Summary of this function goes here
%   Detailed explanation goes here
    trial_length = (ON_dur+(2*OFF_dur))*fs;
    trial = 1;
    
    for i=1:trial_length:length(stimulus)
        stim_clip = stimulus(i:trial_length);
        
        [pks,~] = findpeaks(stim_clip);
        no_of_peaks = length(pks);
        
        stim_frequencies(trial) = no_of_peaks/ON_dur;
        trial=trial+1;
        
    end
end

