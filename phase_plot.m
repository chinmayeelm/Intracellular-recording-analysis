function [gain,phase_lag] = phase_plot(stimulus, response, OFF_dur, ON_dur, fs, period)
    
    stim = stimulus(:,OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    resp = response(:,OFF_dur*fs:(OFF_dur+ON_dur)*fs);
    mean_pos = mean(stim,2);
    [rows, cols] = size(resp);
    stim_freq = [];

    for i=1%:rows
        locs = [];
        zc = [];
        for j = 2:cols-1
          if (stim(i,j-1)<=mean_pos(i) && stim(i,j)>=mean_pos(i) && stim(i,j+1)>mean_pos(i))
              zc(j) = 1;
          end
        end
        [~,locs ]= find(zc==1); 

%         A1 =  subplot(2,1,1); plot(stim(i,:));hold on; plot(locs, stim(i,locs), 'rx'); yline(mean_pos(i)); hold off;

        stim_clips = [];
        resp_clips = [];
        freq = 1 ./ ((diff(locs) / fs));
    %     figure(); plot(stim_freq);
        if mod(length(freq),2)==0
            L=length(freq)-1;
        else
            L=length(freq);
        end
        
        stim_freq(i,:) = sgolayfilt(freq, 1,L);
    %     figure(); plot(stim_freq);
        for k= 2:length(locs)
             stim_clips = stim(i,locs(k-1):locs(k)); %figure(); subplot(2,1,1); plot(stim_clips) ;
             resp_clips = resp(i, locs(k-1):locs(k));% subplot(2,1,2); plot(resp_clips);

             %stim_freq(k-1) = fft_stim(stim_clips, fs, length(stim_clips));
             max_FR(i, k-1) = max(resp_clips);

             [r(i,k-1,:),lags(i,k-1,:)] = xcorr(stim_clips(k-1,:)- mean(stim_clips(k-1,:)), resp_clips(k-1,:)-mean(resp_clips(k-1,:)), 'coeff');
            [~, ind] = max(r(i,k-1,:));


            time_lag = abs(lags(i,k-1,ind))/fs;
            phase_lag = (time_lag/period)*360; %in degrees

            gain = (max(resp_clips)-min(resp_clips))/(max(stim_clips)-min(stim_clips));



%             figure();
            polarscatter(phase_lag, gain);hold on;

        end

        
%         title('Useless plot');
end