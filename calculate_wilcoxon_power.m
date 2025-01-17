function power = calculate_wilcoxon_power(alpha, num_samples, num_simulations, effect_size)
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

        % Perform Wilcoxon rank-sum test (alternative to Nemenyi for two groups)
        p = signrank(group1, group2);

        % Check if the p-value is less than the significance level
        if p < alpha
            rejections = rejections + 1;
        end
    end

    % Estimate power
    power = rejections / num_simulations;
end
