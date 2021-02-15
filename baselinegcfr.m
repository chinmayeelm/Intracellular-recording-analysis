ON_dur = 10;
OFF_dur = 5; 
fs = 10000;
stim_start = OFF_dur*fs;
stim_stop = (ON_dur+OFF_dur)*fs;

% idx = 3;

% stim_all = cell2mat(M1N3_sqr{:, 'antennal_movement'});
% stim = M1N3_sqr.antennal_movement{idx}(:, start:stop);

% resp_all = cell2mat(M1N3_sqr{:, 'gcfr'});
% 
% before_stim = resp_all(:, 1:stim_start+1);
% after_stim =  resp_all(:, stim_stop+1:end);
% 
% before_after = [before_stim after_stim];
% imagesc(before_after);
% spectrogram(before_after);

[rows, ~] = size(M1N3_sqr);
time = M1N3_sqr.time{rows};
for table_row = 1:rows
    
    
    resp = M1N3_sqr.gcfr{table_row};
    [trials,~] = size(resp);
    resp_before = resp(:,1:stim_start+1);
    resp_after = resp(:, stim_stop+1:end);
    
    figure;
    trials= 5;
    for i=1:trials
        subplot(trials,1,i); plot(time(1:length(resp_before)),resp_before(i,:)); hold on;
        plot(time(1:length(resp_before)),resp_after(i,:));
        
        ylabel 'GCFR';
        xlabel 'time (s)';
        title (strcat(M1N3_sqr.stim_name{table_row}, 'Hz', "-", num2str(i)));
        
        
        
    end

end