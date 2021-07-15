function reshaped_data = reshape_data(data, single_trial_length, no_of_protocols, no_of_trials)
%     reshaped_data = (reshape(data, [single_trial_length, no_of_protocols*no_of_trials]))';
    reshaped_data = (reshape(data, [single_trial_length, no_of_trials]))';
    
end