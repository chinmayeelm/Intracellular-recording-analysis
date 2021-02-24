ON_dur = 10;
OFF_dur = 3; 
fs = 10000;
% idx = 17;
start = fs*OFF_dur; stop = (ON_dur+OFF_dur)*fs;

for z=1:8
        
        idx = z;
%         time = T_sqr_27_01.time(idx);
        period = (round(T_sqr_27_01.stim_period(idx),4))*fs;
%         period = T_sqr_27_01.stim_period(idx)*fs;
        t = 0:period;

        stim = T_sqr_27_01.antennal_movement{idx}(1, start:stop);
        resp = T_sqr_27_01.avg_gcfr(idx, start:stop);

        stim = stim - mean(stim,2);

         mean_pos = mean(stim, 2);
            [rows, cols] = size(resp);

            for i=1
                locs = [];
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
                stim_clips = [];
                resp_clips = [];

%                 figure;
%                 within_trial_cycle = 1;
%                 for k= locs(1)%:period:length(stim)-period
                     k = locs(1);
                     disp(k+period);
                     
                     if k+period > ON_dur*fs
                         figure;
                         t = linspace(0,10,ON_dur*fs+1);
                         A1 = subplot(2,1,1); plot(t, resp(1,:)/ max(resp(1,:))); hold on;
%                          A1 = subplot(2,1,1); plot(resp); hold on;
                         ylabel 'GCFR';
                         title ('One cycle across different stimuli');

                         A2 = subplot(2,1,2); plot(t, stim(1,:)); hold on;
%                          A2 = subplot(2,1,2); plot(stim);
                         ylabel 'Antennal movement'
                         xlabel 'time (s)'
                         
                         linkaxes([A1,A2],'x');
                         
                         continue
                     end
                         
                     stim_clips  = stim(i, k:k+period);
                     resp_clips = resp(i, k:k+period);
                     figure;
                     t = linspace(0,period/fs,length(resp_clips));
                     A1 = subplot(2,1,1); plot(t, resp_clips(1,:)/ max(resp_clips(1,:))); hold on;
                     ylabel 'GCFR';
                     title ('One cycle across different stimuli');

                     A2 = subplot(2,1,2); plot(t, stim_clips(1,:)); hold on;
                     ylabel 'Antennal movement'
                     xlabel 'time (s)'
        %              fig_name = sprintf('T_sqr_27_01_%0.0001fs',T_sqr_27_01.stim_period(idx))
%                      within_trial_cycle = within_trial_cycle +1;
                     
                     linkaxes([A1,A2],'x');
%                      xlim([0 4]);
                     
                     

%                 end 
        %         figure;
        %         subplot(2,1,1); stdshade(resp_clips, 0.6, [0.4660 0.6740 0.1880], t); hold on;
        %          ylabel 'GCFR';
        % %          
        %         subplot(2,1,2); stdshade(stim_clips, 0.6, [0.6, 0.2,0], t); hold on;
        %          ylabel 'Antennal movement'
        %          xlabel 'time (s)'

        %          filename = sprintf('T_sqr_27_01_%0.0001fs_t%d.png', T_sqr_27_01.stim_period(idx), i);
        %          saveas(gcf, filename, 'png');
            end
    
end
                 

