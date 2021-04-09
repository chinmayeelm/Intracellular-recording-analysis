ON_dur = 10;
OFF_dur = 3;
fs = 10000;
% idx = 9;
[t_rows,~] = size(T_sqr_M2N4T2);
for idx=1:4
    start = fs*OFF_dur; stop = (ON_dur+OFF_dur)*fs;
    % time = T_sqr_M2N4T2.time(idx);
    % t = time{1,1}(1:T_sqr_M2N4T2.single_trial_length(5));
    period = (round(T_sqr_M2N4T2.stim_period(idx),4))*fs;
    %     t = 0:period;
    % trial_num = 1;
    
    stim = T_sqr_M2N4T2.antennal_movement{idx}(:, start:stop);
    % stim(4:5,:) = [];
    resp = T_sqr_M2N4T2.gcfr{idx}(:, start:stop);
    % resp(4:5,:) = [];
    
    stim = stim - mean(stim,2);
    % amp = max()
    
    mean_pos = mean(stim, 2);
    [rows, cols] = size(resp);
    %     figure;
    
    
    mean_resp_clips = [];
    mean_stim_clips = [];
    
    for i=1:rows
        locs = [];
        zc = [];
        for j = 2:cols-1
            if (stim(i,j-1)>= 0 && stim(i,j)<= 0 && stim(i,j+1)< 0)
                % if (stim(i,j-1)>=mean_pos(i) && stim(i,j)<=mean_pos(i) )%&& stim(i,j+1)>mean_pos(i))
                %             if (stim(i,j)<=mean_pos(i) && stim(i,j+1)>mean_pos(i))
                zc(j) = 1;
            end
        end
        
        [~,locs ]= find(zc==1);
        %         figure; plot(stim(i,:)); hold on; plot(locs, stim(i,locs), 'rx');
        %         figure; plot(resp(i,:));
        stim_clips = [];
        resp_clips = [];
        
        %         figure;
        within_trial_cycle = 1;
        
        for k= locs(1):period:length(stim)-period
            if k+period<ON_dur*fs
                stim_clips(within_trial_cycle,:)  = stim(i, k:k+period);
                resp_clips(within_trial_cycle,:) = resp(i, k:k+period);
                %                 plot(resp_clips); hold on;
                within_trial_cycle = within_trial_cycle +1;
            end
            
        end
        
        %          t = linspace(0,period/fs,period+1);
        %          a1 = subplot(2,1,1); plot(t, mean(resp_clips)); hold on;
        %          ylabel 'Avg GCFR';
        %          title (strcat("stim period=", num2str(T_sqr_M2N4T2.stim_period(idx)), " s"));
        %
        %
        %          a2 = subplot(2,1,2); plot(t, mean(stim_clips)); hold on;
        %          ylabel 'Antennal movement'
        %          xlabel 'time (s)'
        %
        %          linkaxes([a1,a2], 'x');
        
        mean_resp_clips(i,:) = mean(resp_clips,1);
        mean_stim_clips(i,:) = mean(stim_clips,1);
    end
    
    if sum(isnan(mean_stim_clips),'all')==0 && ~isempty(mean_resp_clips)
        figure;
        t = linspace(0,period/fs,period+1);
        a1 = subplot(2,1,1); stdshade(mean_resp_clips, 0.6, [0.4660 0.6740 0.1880], t); hold on;
        ylabel 'Avg. GCFR';
        title (strcat("Avg. of cycles within trials.   ", "stim period=", num2str(T_sqr_M2N4T2.stim_period(idx)), " s   ", T_sqr_M2N4T2.filename(idx), "_", T_sqr_M2N4T2.date(idx)), 'interpreter', 'none');
        
        %
        a2 = subplot(2,1,2); stdshade(mean_stim_clips, 0.6, [0.6, 0.2,0], t); hold on;
        ylabel 'Antennal movement'
        xlabel 'time (s)'
        
        linkaxes([a1,a2], 'x');
        filename = strcat("Avg of all cycles_", "stim period", join(split(num2str(T_sqr_M2N4T2.stim_period(idx)),'.')), " s_", T_sqr_M2N4T2.filename(idx),"t",num2str(i));
        saveas(gcf, filename, 'png');
    else
        figure;
        t = linspace(0,ON_dur,ON_dur*fs+1);
        a1 = subplot(2,1,1); stdshade(resp, 0.6, [0.4660 0.6740 0.1880], t); hold on;
        ylabel 'Avg. GCFR';
        title (strcat("Avg. of cycles within trials.   ", "stim period=", num2str(T_sqr_M2N4T2.stim_period(idx)), " s   ", T_sqr_M2N4T2.filename(idx), "_", T_sqr_M2N4T2.date(idx)), 'interpreter', 'none');
        
        %
        a2 = subplot(2,1,2); stdshade(stim, 0.6, [0.6, 0.2,0], t); hold on;
        ylabel 'Antennal movement'
        xlabel 'time (s)'
        
        linkaxes([a1,a2], 'x');
        filename = strcat("Avg of all cycles_", "stim period", join(split(num2str(T_sqr_M2N4T2.stim_period(idx)),'.')), " s_", T_sqr_M2N4T2.filename(idx),"t",num2str(i));
        saveas(gcf, filename, 'png');
    end
    
end


