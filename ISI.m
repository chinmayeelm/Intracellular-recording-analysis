ON_dur = 10;
OFF_dur = 3; 
fs = 10000;
idx = 1;

start = fs*OFF_dur; stop = (ON_dur+OFF_dur)*fs;
period = (round(M1_N6_T2.stim_period(idx),4))*fs;

stim = M1_N6_T2.antennal_movement{idx}(:, start:stop);
resp = M1_N6_T2.raster{idx}(:, start:stop);
gcfr = M1_N6_T2.gcfr{idx}(:,start:stop);

stim = stim - mean(stim,2);

mean_pos = mean(stim, 2);
[rows, cols] = size(resp);

for i=1:rows
    locs = [];
    zc = [];
    for j = 2:cols-1
      if (stim(i,j-1)>= 0 && stim(i,j)<= 0 && stim(i,j+1)< 0)
          zc(j) = 1;
      end
    end

    [~,locs ]= find(zc==1);
    stim_clips = [];
    resp_clips = [];
    gcfr_clips = [];

%         figure;
    within_trial_cycle = 1;
    for k= locs(1) %:period :length(stim)-period
         stim_clips(within_trial_cycle,:)  = stim(i, k:k+period);
         resp_clips(within_trial_cycle,:) = resp(i, k:k+period);
         gcfr_clips(within_trial_cycle, :) = gcfr(i, k:k+period);
         
         spike_locs = find(resp_clips(within_trial_cycle,:));
         isi = diff(spike_locs);

%              figure;
             
             subplot(3,1,1); plot(spike_locs(2:end)/fs, isi); hold on; 
             ylabel 'ISI';
             title (strcat("stim period = ", num2str(M1_N6_T2.stim_period(idx)), " s"));
             xlim ([0 period/fs]);
             
             t = linspace(0,period/fs,period+1);
             subplot(3,1,2); plot(t, gcfr_clips(within_trial_cycle, :)); hold on;
             ylabel 'GCFR';
              
                 
             subplot(3,1,3); plot(t, stim_clips(within_trial_cycle,:)); hold on;
             ylabel 'Antennal movement'
             xlabel 'time (s)'
%              fig_name = sprintf('M1_N6_T2_%0.0001fs',M1_N6_T2.stim_period(idx))
         within_trial_cycle = within_trial_cycle +1;

    end 

end
