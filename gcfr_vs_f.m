function [gcfr_at_f, f] = gcfr_vs_f(stimulus, response, sampling_rate, ON_duration, OFF_duration)
start_point = OFF_duration * sampling_rate;
end_point = (OFF_duration + ON_duration) * sampling_rate;
stimulus = stimulus(start_point:end_point);
stimulus = stimulus - mean(stimulus);
response = response(start_point:end_point);
crossing_level = 0.05 * (max(stimulus) - min(stimulus));

levels_crossed_at = find(diff(stimulus > -crossing_level & stimulus < crossing_level));
f = 1 ./ (diff(levels_crossed_at) / sampling_rate);
gcfr_at_f = zeros(1,length(f));
for i = 1:length(levels_crossed_at)-1
    gcfr_at_f(i) = max(response(levels_crossed_at(i):levels_crossed_at(i+1)));
end
end