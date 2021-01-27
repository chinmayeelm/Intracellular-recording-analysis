ON_dur = 10;
OFF_dur = 3; 
fs = 10000;
stim_start = OFF_dur*fs;
stim_stop = (ON_dur+OFF_dur)*fs;

% idx = 3;

% stim_all = cell2mat(T_sqr{:, 'antennal_movement'});
% stim = T_sqr.antennal_movement{idx}(:, start:stop);

% resp_all = cell2mat(T_sqr{:, 'gcfr'});
% 
% before_stim = resp_all(:, 1:stim_start+1);
% after_stim =  resp_all(:, stim_stop+1:end);
% 
% before_after = [before_stim after_stim];
% imagesc(before_after);
% spectrogram(before_after);
 time = T_sqr.time{idx};

[idx, ~] = size(T_sqr);
for table_row = 1:idx
    
    
    resp = T_sqr.gcfr{idx};
    [trials,~] = size(resp);
    resp_before = resp(:,1:stim_start+1);
    resp_after = resp(:, stim_stop+1:end);
    
    figure;
    for i=1:trials
        subplot(trials,1,i); plot(time(1:length(resp_before)),resp_before(i,:)); hold on;
        plot(time(1:length(resp_before)),resp_after(i,:));
        
        ylabel 'GCFR';
        xlabel 'time (s)';
        title (strcat(T_sqr.stim_name{table_row}, 'Hz', "-", num2str(i)));
        
        
        
    end

end