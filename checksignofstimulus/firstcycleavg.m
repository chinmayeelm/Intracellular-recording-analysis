ON_dur = 10;
OFF_dur = 3;
fs = 10000;
% idx = 17;
start = fs*OFF_dur; stop = (ON_dur+OFF_dur)*fs;

for idx=1:9
    
    
    %         time = T_sqr.time(idx);
    period = (round(T_sqr.stim_period(idx),4))*fs;
    %         period = T_sqr.stim_period(idx)*fs;
    %     t = 0:period;
    
    stim = T_sqr.antennal_movement{idx}(:, start:stop);
    resp = T_sqr.gcfr{idx}(:, start:stop);
    
    stim = stim - mean(stim,2);
    
    mean_pos = mean(stim(1,:), 2);
    [rows, cols] = size(resp);
    
    stim_clips = [];
    resp_clips = [];
%     figure;
    for i=1
        
        zc = [];
        for j = 2:cols-1
            if (stim(i,j-1)>= 0 && stim(i,j)<= 0 && stim(i,j+1)< 0)
                % if (stim(i,j-1)>=mean_pos(i) && stim(i,j)<=mean_pos(i) )%&& stim(i,j+1)>mean_pos(i))
                %             if (stim(i,j)<=mean_pos(i) && stim(i,j+1)>mean_pos(i))
                zc(j) = 1;
            end
        end
        
        [~,locs]= find(zc==1);
        %         figure; plot(stim(i,:)); hold on; plot(locs, stim(i,locs), 'rx');
        %         figure; plot(resp(i,:));
        
        
        %                 figure;
        %                 within_trial_cycle = 1;
        %                 for k= locs(1)%:period:length(stim)-period
        k = locs(1);
        disp(k+period);
        
        
        if k+period > ON_dur*fs
            %             figure;
            t = linspace(0,10,ON_dur*fs+1);
            
            A1 = subplot(2,1,1); 
            plot(t, mean(resp)./max(mean(resp))); hold on;
            %stdshade(resp, 0.6, [0.4660 0.6740 0.1880], t); hold on;
            ylabel 'Avg. GCFR';
            title (strcat("One stimulus cycle", "stim period=", num2str(T_sqr.stim_period(idx)), " s   ", T_sqr.filename(idx)," ", T_sqr.date(idx)), 'interpreter', 'none');
            
            
            A2 = subplot(2,1,2); 
            plot(t,stim(i,:)); hold on;
%             stdshade(stim, 0.6, [0.6, 0.2,0], t); hold on;
            ylabel 'Antennal movement'
            xlabel 'time (s)'
            
            linkaxes([A1,A2],'x');
            
%             filename = strcat("First cycle avg_", "stim period", join(split(num2str(T_sqr.stim_period(idx)), '.')), " s_", T_sqr.filename(idx),"t",num2str(i));
%             saveas(gcf, filename, 'png');
            
            break
        end
        
        stim_clips  = stim(i, k:k+period);
        resp_clips = mean(resp(:, k:k+period));
        %         figure;
        %         t = linspace(0,period/fs,length(resp_clips));
        %         A1 = subplot(2,1,1); plot(t, resp_clips(1,:)/ max(resp_clips(1,:))); hold on;
        %         ylabel 'GCFR';
        %         title (strcat("One stimulus cycle", "stim period=", num2str(T_sqr.stim_period(idx)), " s   ", T_sqr.filename(idx)), 'interpreter', 'none');
        %
        %         A2 = subplot(2,1,2); plot(t, stim_clips(1,:)); hold on;
        %         ylabel 'Antennal movement'
        %         xlabel 'time (s)'
        %         %              fig_name = sprintf('T_sqr_27_01_%0.0001fs',T_sqr.stim_period(idx))
        %         %                      within_trial_cycle = within_trial_cycle +1;
        %         linkaxes([A1,A2],'x');
        %                      xlim([0 4]);
        
        
        
    end
    
    t = linspace(0,period/fs,length(resp_clips));
    subplot(2,1,1); plot(t, resp_clips./max(resp_clips)); hold on;
%     stdshade(resp_clips, 0.6, [0.4660 0.6740 0.1880], t); hold on;
    ylabel 'GCFR';
%     title (strcat("One stimulus cycle", "stim period=", num2str(T_sqr.stim_period(idx)), " s   ", T_sqr.filename(idx)," ", T_sqr.date(idx)), 'interpreter', 'none');
    %
    subplot(2,1,2); plot(t,stim_clips); hold on;
%     stdshade(stim_clips, 0.6, [0.6, 0.2,0], t); hold on;
    ylabel 'Antennal movement'
    xlabel 'time (s)'
    
%     filename = strcat("First cycle avg_", "stim period", join(split(num2str(T_sqr.stim_period(idx)),'.')), " s_", T_sqr.filename(idx),"t",num2str(i));
%     saveas(gcf, filename, 'png');
    
end




