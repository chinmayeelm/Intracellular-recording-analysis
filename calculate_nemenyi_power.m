function power = calculate_nemenyi_power(alpha, num_samples, num_simulations, effect_size)
    % Parameters
    % alpha: significance level (e.g., 0.05)
    % num_samples: number of samples in each group (e.g., 5)
    % num_simulations: number of simulations to estimate power
    % effect_size: the difference between the means of two groups to detect

    % Initialize counters
    rejections = 0;

    for i = 1:num_simulations
        % Generate data for two groups with given effect size
        group1 = randn(num_samples, 1);               % Control group
        group2 = randn(num_samples, 1) + effect_size; % Test group with effect size

        % Combine data and ranks
        data = [group1; group2];
        ranks = tiedrank(data);
        
        % Assign ranks to each group
        ranks1 = ranks(1:num_samples);
        ranks2 = ranks(num_samples + 1:end);
        
        % Compute mean rank for each group
        mean_rank1 = mean(ranks1);
        mean_rank2 = mean(ranks2);
        
        % Calculate test statistic for Nemenyi (absolute difference of mean ranks)
        q_statistic = abs(mean_rank1 - mean_rank2);
        
        % Critical value for Nemenyi test (for two groups, simplified)
        % Approximated for small groups; adjust with more groups
        critical_value = sqrt((12 * num_samples * (num_samples + 1)) / (2 * num_samples)) * norminv(1 - alpha / 2);
        
        % Check if the test statistic exceeds the critical value
        if q_statistic > critical_value
            rejections = rejections + 1;
        end
    end

    % Estimate power
    power = rejections / num_simulations;
end
