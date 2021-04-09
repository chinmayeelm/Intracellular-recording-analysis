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

[rows, ~] = size(T_sqr);
time = T_sqr.time{rows};
for table_row = 1:rows
    
    
    resp = T_sqr.gcfr{table_row};
    [trials,~] = size(resp);
    resp_before = resp(:,1:stim_start+1);
    resp_after = resp(:, stim_stop+1:end);
    
    figure;
%     trials= 5;
    for i=1:trials
        subplot(trials,1,i); plot(time(1:length(resp_before)),resp_before(i,:)); hold on;
        plot(time(1:length(resp_before)),resp_after(i,:));
        
        ylabel 'GCFR';
        xticklabels([]);
        
        if i==1
            title (strcat("Baseline GCFR before and after stimulus. ", T_sqr.stim_name{table_row}, 'Hz', T_sqr.filename(table_row), "_", T_sqr.date(table_row)), 'Interpreter', 'none');
        end
  
        if i==trials
            xticklabels('auto');
            xlabel 'time (s)';
        end
      
    end
    filename= strcat("Baseline GCFR before and after stimulus_", join(split(T_sqr.stim_name{table_row},'.')),'Hz_', T_sqr.filename{table_row});
    saveas(gcf,filename,'png');

end