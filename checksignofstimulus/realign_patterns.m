function shifted_stim_patterns = realign_patterns(stim_patterns,shift_time,fs)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
for pattern_idx = 1:length(shift_time)
    shift_by_samples = shift_time(pattern_idx)*fs;
    shifted_stim_patterns = circshift(stim_patterns, shift_by_samples);
end

