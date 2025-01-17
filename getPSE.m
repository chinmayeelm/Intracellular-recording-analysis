function stimulus_prior = getPSE(stimulus, stim_window, fs, NstimPriors)
%GETPSE Returns required number of random stimulus patterns of length
%stim_window

pattern_length = stim_window*fs;
stimulus_prior = zeros([NstimPriors pattern_length]);
rr = randi([1 (length(stimulus)-pattern_length)],1,NstimPriors);
parfor j=1:NstimPriors
    stimulus_prior(j,:) = stimulus(randperm(size(stimulus,1),1),rr(j):rr(j)+pattern_length-1);
end

end