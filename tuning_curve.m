function [stim_freq, max_FR]  = tuning_curve(stim, resp, fs, ON_dur, OFF_dur)

    stim = mean(stim(:,OFF_dur*fs:(OFF_dur+ON_dur)*fs));
    resp = resp(:,OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    mean_pos = mean(stim);
    [rows, cols] = size(resp);
    stim_freq = [];

    for i=1:rows
        locs = [];
        zc = [];
        for j = 2:cols-1
          if (stim(j-1)<=mean_pos && stim(j)>=mean_pos && stim(j+1)>mean_pos)
              zc(j) = 1;
          end
        end
        [~,locs ]= find(zc==1); 

%         A1 =  subplot(2,1,1); plot(stim(i,:));hold on; plot(locs, stim(i,locs), 'rx'); yline(mean_pos(i)); hold off;

        stim_clips = [];
        resp_clips = [];
        freq = 1 ./ ((diff(locs) / fs));
    %     figure(); plot(stim_freq);
        stim_freq = sgolayfilt(freq, 1,51);
    %     figure(); plot(stim_freq);
        for k= 2:length(locs)
             stim_clips = stim(locs(k-1):locs(k)); %figure(); subplot(2,1,1); plot(stim_clips) ;
             resp_clips = resp(i, locs(k-1):locs(k));% subplot(2,1,2); plot(resp_clips);

             %stim_freq(k-1) = fft_stim(stim_clips, fs, length(stim_clips));
             max_FR(i, k-1) = max(resp_clips);

    %          [r(k-1,:),lags(k-1,:)] = xcorr(stim_clips(k-1,:)- mean(stim_clips(k-1,:)), resp_clips(k-1,:)-mean(resp_clips(k-1,:)), 'coeff');
        end

        size(stim_freq)
        size(max_FR)
       
    end
%     A2 = subplot(2,1,2); plot(stim_freq, mean(max_FR));
%     figure();
%     [lineOut, fillOut] = stdshade(max_FR,0.2,'k',stim_freq);
%     lineOut.LineWidth = 0.5;
%     ylabel('Normalised firing rate (spike/s)');
%     xlabel('Stimulus frequency (Hz)');
    %         time_lag = abs(lags)./fs;
    %         phase_lag = (time_lag/diff(locs)).*360; %in degrees
    % 
    %         gain = (max(stim_clips)-min(stim_clips))/(max(resp_clips)-min(resp_clips));
    %         
    %         
    % 
    %         figure();
    %         polarscatter(phase_lag, gain);
    %         title('Useless plot');

    %     linkaxes([A1, A2],'x');
 
end
