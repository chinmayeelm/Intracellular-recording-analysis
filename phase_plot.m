function [corr_coef_val, vs, phase] = phase_plot(stimulus, response, raster,start_stim, stop_stim, fs, period)

stimulus = stimulus(:,start_stim : stop_stim);
resp = response(:,start_stim : stop_stim);
raster = raster(:,start_stim : stop_stim);

zero_ind = find(raster==0);
raster(zero_ind) = NaN;

stim = stimulus - mean(stimulus,2);

[ntrials,nsamples] = size(stim);
mean_positions = mean(stim, 2);
resp = resp - mean(resp,2);
% figure;
%     plot(stim(m,:)); hold on;
%     plot(raster(m,:).*4e-3, 'r|'); hold off;
phase = [];

figure;
for trial = 1:ntrials
    
    cos_theta = [];
    sin_theta = [];
    count = 0;
    mean_pos  =mean_positions(trial);
    
    locs = [];
    zc = [];
    for sample = 2:nsamples-1
        if (stim(trial,sample-1)<=mean_pos && stim(trial,sample)>=mean_pos && stim(trial,sample+1)>mean_pos)
            zc(sample) = 1;
        end
    end
    [val,locs ]= find(zc==1);
    
    %     figure; plot(stim(i,:));hold on; plot(locs, stim(i,locs), 'rx'); yline(mean_pos); hold off;
    
    for k= 2:length(locs)
        
        
        stim_clips = stim(trial, locs(k-1):locs(k)); %figure(); subplot(2,1,1); plot(stim_clips);
        resp_clips = resp(trial, locs(k-1):locs(k)); %subplot(2,1,2); plot(resp_clips);
        raster_clips = raster(trial, locs(k-1):locs(k));
        
        %         subplot(2,1,2); plot(raster_clips, 0, '|'  );
        
        if period > 0.2
            
            %             count = count + 1;
            %             coef = corrcoef(stim_clips, resp_clips);
            %             corr_coef(count)  = coef(1,2);
            
            
            
        elseif period <=0.2
            %             figure; plot(stim_clips);
            %             figure; plot(raster_clips);
            %             [~,stim_peak_ind] = max(stim_clips);
            ind = find(raster_clips == 1);
            if length(ind)>=1
                count = count+1;
                theta = (2*pi*ind(1)/fs)/period;
                phase(count) = theta; %get_spike_phase(ind(1), fs, period);
                
                
                cos_theta(count) = cos(theta);
                sin_theta(count) = sin(theta);
                
                coef = corrcoef(stim_clips, resp_clips);
                corr_coef(count)  = coef(1,2);
            end
            
        end
    end
    
    if period > 0.2
        vs = NaN;
        phase = NaN;
        %         corr_coef_val(i) = mean(corr_coef);
        coef = corrcoef(stim(trial,:), resp(trial,:));
        corr_coef_val(trial) = coef(1,2);
        plot(period, corr_coef_val(trial), '.'); hold on;
        
    elseif period <=0.2
        sum(cos_theta);
        sum(sin_theta);
        vs(trial) = (sqrt(((sum(cos_theta))^2)+((sum(sin_theta))^2)))/count;
        p(trial) = mean(phase);
        polarscatter(p(trial),vs(trial), 100, '.'); hold on;
        rlim([0 1])
        title("frequency = " + 1/period + "Hz");
        corr_coef_val(trial) = mean(corr_coef);
    end
    
end
end
