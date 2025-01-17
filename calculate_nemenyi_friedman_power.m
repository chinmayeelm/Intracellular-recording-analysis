function power = calculate_nemenyi_friedman_power(alpha, num_groups, num_samples, num_simulations, effect_size)
    % Parameters
    % alpha: significance level (e.g., 0.05)
    % num_groups: number of groups (e.g., 6)
    % num_samples: number of samples per group (e.g., 5)
    % num_simulations: number of simulations to estimate power
    % effect_size: the difference in mean rank between pairs of groups

    % Initialize counters for significant pairwise comparisons
    rejections = 0;
    total_comparisons = 0;

    for i = 1:num_simulations
        % Generate data with an effect size applied to one of the group pairs
        data = randn(num_samples, num_groups);  % Random normal data

        % Introduce an effect size to simulate difference (e.g., between group 1 and 2)
        data(:, 1) = data(:, 1) + effect_size;  % Shift group 1 up by effect size
        
        % Perform Friedman test
        [p_friedman, ~, stats] = friedman(data, 1, 'off'); % 'off' suppresses the display

        % If Friedman test is significant, proceed with Nemenyi post-hoc test
        if p_friedman < alpha
            % Compute ranks and mean ranks for each group
            ranks = tiedrank(data')';
            mean_ranks = mean(ranks, 1);

            % Perform pairwise Nemenyi test comparisons
            for j = 1:num_groups - 1
                for k = j + 1:num_groups
                    % Compute the q-statistic for pairwise comparisons
                    q_statistic = abs(mean_ranks(j) - mean_ranks(k));
                    
                    % Calculate critical value for Nemenyi test
                    critical_value = sqrt((num_groups * (num_groups + 1)) / (6 * num_samples)) * ...
                                     norminv(1 - alpha / (2 * num_groups * (num_groups - 1)));

                    % Check if q-statistic exceeds critical value
                    if q_statistic > critical_value
                        rejections = rejections + 1;
                    end
                    total_comparisons = total_comparisons + 1;
                end
            end
        end
    end

    % Estimate power as the proportion of significant pairwise comparisons
    power = rejections / total_comparisons;
end
