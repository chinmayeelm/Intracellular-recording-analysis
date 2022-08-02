function dejittered_STA = dejittering(stimulus,raster, window, fs)
%Implementation of spike dejittering as per Aldworth et al.
%   Detailed explanation goes here

stim_patterns = extract_stim_patterns(stimulus,raster, window, fs);

%     STA = mean(stim_patterns,1);


sigma = 3e-3; %3 ms
lmin  = 1e-4; % 0.1 ms
step = 1e-4;


err = 1;
while (err <= 1e-6)
    shift_vals = (lmin:step:sigma);
    
    stim_patterns_cov_mat = cov(stim_patterns);
    
    mean_waveform = mean(stim_patterns,1);
    plot(mean_waveform); title(err);
    
    shift_time = [];
    
    
    for stim_pattern_idx = 1:size(stim_pattterns,1)
        
        distance = zeros(1, length(shift_vals));
        for shift_iter = 1:length(shift_vals)
            
            stim_pattern_shifted = circshift(stim_patterns, shift_vals(shift_iter));
            t = shift_vals(shift_iter);
            
            residual_waveform = stim_pattern_shifted - mean_waveform;
            distance(shift_iter) = 0.5*((residual_waveform*stim_patterns_cov_mat*residual_waveform') + (t/sigma)^2);
            
        end
        
        [~,idx] = min(distance);
        shift_time(stim_pattern_idx) = shift_vals(idx);
        
    end
    
    stim_patterns_prev = stim_patterns;
    stim_patterns = realign_patterns(stim_patterns, shift_time, fs);
    sigma = var(shift_time);
    
    err = (mean(var(stim_patterns_prev,0,1))-mean(var(stim_patterns,0,1)))/mean(var(stim_patterns_prev,0,1));
end

dejittered_STA = mean(stim_patterns);

end



