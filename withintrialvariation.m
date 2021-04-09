ON_dur = 10;
OFF_dur = 3;
fs = 10000;
idx = 1;
start = fs*OFF_dur; stop = (ON_dur+OFF_dur)*fs;
% time = T_sqr.time(idx);
% t = time{1,1}(1:T_sqr.single_trial_length(5));

for idx=1:18
    period = (round(T_sqr.stim_period(idx),4))*fs;
    t = 0:period;
    % trial_num = 1;
    
    stim = T_sqr.antennal_movement{idx}(:, start:stop);
    % stim(4:5,:) = [];
    resp = T_sqr.gcfr{idx}(:, start:stop);
    % resp(4:5,:) = [];
    
    stim = stim - mean(stim,2);
    % amp = max()
    
    mean_pos = mean(stim, 2);
    [rows, cols] = size(resp);
    %     figure;
    
    %%
    
    
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
        
        k = locs(1);
        
        if k+period > ON_dur*fs
            figure;
            t = linspace(0,ON_dur,ON_dur*fs+1);
            % A1 = subplot(2,1,1); plot(t, resp(:,:)/ max(resp(:,:))); hold on;
            A1 = subplot(2,1,1); plot(t, resp); hold on;
            ylabel 'GCFR';
            title (strcat("Within trial variation.      ", "stim period=", num2str(T_sqr.stim_period(idx)), " s   ", T_sqr.filename(idx)," ", T_sqr.date(idx)), 'interpreter', 'none');
            
            A2 = subplot(2,1,2); plot(t, stim); hold on;
            %                          A2 = subplot(2,1,2); plot(stim);
            ylabel 'Antennal movement'
            xlabel 'time (s)'
            
            linkaxes([A1,A2],'x');
            
            filename = strcat("Within trial variation_", "stim period", join(split(num2str(T_sqr.stim_period(idx)),'.')), " s_", T_sqr.filename(idx),"t",num2str(i));
            saveas(gcf, filename, 'png');
            
            break
        end
        
        figure;
        within_trial_cycle = 1;
        for k= locs(1):period:length(stim)-period
            stim_clips(within_trial_cycle,:)  = stim(i, k:k+period);
            resp_clips(within_trial_cycle,:) = resp(i, k:k+period);
            
            t = linspace(0,period/fs,period+1);
            %         figure;
            subplot(2,1,1); plot(t, resp_clips(within_trial_cycle,:)); hold on;
            ylabel 'GCFR';
            title (strcat("Within trial variation.      ", "stim period=", num2str(T_sqr.stim_period(idx)), " s   ", T_sqr.filename(idx)," ", T_sqr.date(idx)), 'interpreter', 'none');
            
            
            subplot(2,1,2); plot(t, stim_clips(within_trial_cycle,:)); hold on;
            ylabel 'Antennal movement'
            xlabel 'time (s)'
            
            within_trial_cycle = within_trial_cycle +1;
            
        end
        %
        %     figure;
        %     subplot(2,1,1); stdshade(resp_clips, 0.6, [0.4660 0.6740 0.1880], t); hold on;
        %     ylabel 'GCFR';
        %     %
        %     subplot(2,1,2); stdshade(stim_clips, 0.6, [0.6, 0.2,0], t); hold on;
        %     ylabel 'Antennal movement'
        %     xlabel 'time (s)'
        
        filename = strcat("Within trial variation_", "stim period", join(split(num2str(T_sqr.stim_period(idx)),'.')), " s_", T_sqr.filename(idx),"t",num2str(i));
        saveas(gcf, filename, 'png');
    end
    
    
end

