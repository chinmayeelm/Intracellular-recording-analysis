function [jitter] = jitter_SD(raster, stim_window, spike_window, min_no_trials_required, fs)
%JITTER_SD Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(raster);

locs_ref = find(raster(m,:)==1);

% ref = locs_ref(1);
jitter = [];

for j = 1:length(locs_ref)

    %     if locs_ref(j) > ref+spike_window
    ref = locs_ref(j);

    if ref>stim_window && (ref+spike_window)< n


        rows_with_multiple_ones = find(sum(raster(:,ref-spike_window : ref+spike_window), 2) > 1);
        if ~isempty(rows_with_multiple_ones)
            disp(rows_with_multiple_ones);
        end
        [rows, locs] = find(raster(:,ref-spike_window : ref+spike_window)==1);
        mult_spike_idx = ismember(rows, rows_with_multiple_ones);
        locs(mult_spike_idx) = [];
        locs_trials = ref+locs-spike_window-1; % 11 to get the actual location from raster


        if length(locs_trials)<min_no_trials_required || std(locs_trials)==0
            % disp(length(locs_trials));
            % disp(std(locs_trials));
            continue;
        end

        SD = (std(locs));
        jitter = [jitter; SD*1000/fs]; % * 1000 to convert to ms
    end
end
end

