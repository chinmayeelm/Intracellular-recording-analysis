function [positions,freq] = getFirstSpikeLocs(stimulus, antennal_position, raster, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
zero_ind = raster==0;
raster(zero_ind) = NaN;

[ntrials,nsamples] = size(stimulus);
mean_pos = mean(stimulus, "all");
positions = [];
freq = [];

zc = [];
meanStim = mean(stimulus,1);
for sample = 2:nsamples-5
    if (meanStim(sample-1)<=mean_pos && meanStim(sample)>=mean_pos && meanStim(sample+5)>mean_pos)
        zc(sample) = 1;
    end
end



[val,locs ]= find(zc==1);

figure; plot(meanStim);hold on; plot(locs, meanStim(locs), 'rx'); yline(mean_pos); hold off;
% figure;
for trial = 1:ntrials




    for k= 2:length(locs)


        stim_clips = antennal_position(trial, locs(k-1):locs(k)); %figure(); subplot(2,1,1); plot(stim_clips);
        
        raster_clips = raster(trial, locs(k-1):locs(k));

        %         subplot(2,1,2); plot(raster_clips, 0, '|'  );
        if length(stim_clips) > 130
            % figure; plot(stim_clips); 
            % yyaxis right; plot(raster_clips, 'r|');
            % pause;
            %             [~,stim_peak_ind] = max(stim_clips);
            ind = find(raster_clips == 1, 1);
            if length(ind)>=1
                positions = [positions; stim_clips(ind)];
                freq = [freq; 1*fs/length(stim_clips)];

            end

        end


    end
end
end
