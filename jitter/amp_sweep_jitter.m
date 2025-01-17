
cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
filepath = readlines('ampsweeplist.txt');
dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);

cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " sweep");
c = jet(nfiles);

high_jitter = [16 14 12 6 1];
all_files = 1:nfiles;
low_jitter = setdiff(all_files,high_jitter);
T_amp_sweep_jitter = table();
amp = [];
jitter = [];

figure; 
for irow = 1:nfiles
    irow
    P = getStructP(dataDirectory(irow), filename(irow),[nan nan],1);
    fs = P(1).fs;

    startPt = P(1).OFF_dur*P(end).fs+1;
    stopPt = (P(1).OFF_dur + P(end).ON_dur)*P(end).fs;
    resp = P(1).raster(:,startPt:stopPt);
    stim = P(1).mean_movement(startPt:stopPt);
    mean_pos = mean(stim(1,:));
    time = linspace(0,P(end).ON_dur,P(end).ON_dur*fs);


    [pks,locs] = findpeaks(stim,"MinPeakDistance",fs/(2*P(1).amp_sweep_frq), "MinPeakHeight",0.01);

    idx = 1;
    while ~isempty(idx)
        idx = [];
        period_at_locs = 1./P(1).amp_sweep_frq;
        diff_locs = diff(locs)/fs;
        for i=1:length(diff_locs)
            if diff_locs(i)<period_at_locs/2
                idx = [idx;i];
            end
        end
        locs(idx+1)=[];
    end



    % figure;
    % plot(time, stim(1,:)); hold on; plot(locs/fs, stim(1, locs), 'rx'); yline(mean_pos); hold off;
    [rr,cc] = find(resp);
    % yyaxis right; plot(cc/fs,rr,'k|');

    stim_clips = [];
    resp_clips = [];
    I_spike_phase = [];
    amp_val = [];
    spike_time_sd = [];

    amp_val_list = [];
    spike_time_list = [];

    I = 1;

       
    % figure;
    % plot(time, stim(1,:)); hold on; plot(locs/fs, stim(1, locs), 'rx'); yline(mean_pos); hold off;
    % yyaxis right; plot(cc/fs,rr,'k|');hold on;


    for k= 2:length(locs)
        stim_clips = stim(1, locs(k-1):locs(k));
        resp_clips = resp(:, locs(k-1):locs(k));


        % yyaxis left; plot((locs(k-1):locs(k))./fs,stim_clips, 'k-');
        

        [row_ind,col_ind] = find(resp_clips == 1);
        % yyaxis right; plot((locs(k-1)+col_ind)./fs, row_ind, 'm|');
        % pause;

        [C, ia, ic] = unique(row_ind, 'stable');
        ind = col_ind(ia);




        if length(col_ind)>=1
            I_spike_time = ind/fs;
            if length(ind) == P(1).complete_trials
                spike_time_sd(I) = (std(I_spike_time))*1000;
                %             I_spike(I) = freq(ind(1)+locs(k-1));%stim_freq(k-1);
                amp_val(I) = max(stim_clips) - min(stim_clips);
                mean_spike_time(I) = mean(I_spike_time)*1000;

                if amp_val(I) > 0.01
                    amp_val_list = [amp_val_list amp_val(I).*ones(1,length(ind))];
                    spike_time_list = [spike_time_list (ind.*1000/fs)'];

                    %                 scatter(amp_val(I).*ones(1,length(ind)), ((ind./fs)'));
                    %                 hold on;
                end

                I=I+1;

            end
        end



    end


    amp = [amp; amp_val'];
    jitter = [jitter; spike_time_sd'];
    T_amp_sweep_jitter.amp(irow) = {amp};
    T_amp_sweep_jitter.jitter(irow) = {jitter};
    % figure;
    scatter(amp_val, spike_time_sd, 10, "o", "MarkerFaceColor",c(irow,:), "MarkerEdgeColor","none","MarkerFaceAlpha",0.8);
    hold on;
    % title(irow)
    % pause;
end

ylabel('Jitter (ms)');
xlabel('Amplitude ({\circ})');
axis padded;
yticks(0:10:50);

%% Spearman's rho

for irow=1:nfiles

    amp = T_amp_sweep_jitter.amp{irow,1};
    jitter = T_amp_sweep_jitter.jitter{irow,1};
    rho = corr(amp,jitter, 'type','Pearson' );

    T_amp_sweep_jitter.rho(irow) = rho;
end
