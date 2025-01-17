cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
filepath = readlines('chirpFullList.txt');
dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);

cellID = replace(join([dataDirectory filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " chirp");
c = parula(nfiles);

T_chirp = table();

for irow = 9%:nfiles

    P = getStructP(dataDirectory(irow), filename(irow),[nan nan],1);

    startPt = P(end).OFF_dur*P(end).fs+1;
    stopPt = (P(end).OFF_dur + P(end).ON_dur)*P(end).fs;

    resp = P(end).raster(:,startPt:stopPt);
    stim = P(end).antennal_movement(1,startPt:stopPt);
    fs = P(end).fs;
    % order = 8;
    % [x,y] = butter(order, P(end).max_chirp_frq*2/(fs/2), 'low');
    % stim = filtfilt(x,y,stim);
    
    stim_ifb = P(end).stim_ifb(1,startPt:stopPt);
    corrcoef(stim,stim_ifb)
    mean_pos = mean(stim);
    f_fit = miscFuncs.getFreqFit(stim_ifb, fs, P(end).ON_dur);
    % [f,t] = instfreq(stim,fs, 'FrequencyLimits',[1 180]);
    % f_fit = fit(t,f,'poly3');
    time = linspace(0,P(end).ON_dur,length(stim));
    freq = f_fit(time); %P(end).inc_frq_chirp_f;
    % figure; plot(f_fit,time,freq);
    T_chirp.freq(irow) = {freq};
    T_chirp.gcfr(irow) = {P(end).gcfr(:,startPt:stopPt)};
%{
    % [c, midlev] = midcross(stim, fs);
    [rr,cc] = find(resp);
    % figure;
    % plot(time, stim(1,:)); hold on; %plot(locs/fs, stim(locs), 'rx'); yline(mean_pos);
    % plot(c, midlev, 'rx');
    % yyaxis right; plot(cc/fs,rr,'k|')
    % hold off;

    stim_smoothened = sgolayfilt(stim_ifb,1,51);

    [pks,locs] = findpeaks(stim_smoothened,"MinPeakDistance",fs/(2*max(freq)), "MinPeakHeight",0.05);

    idx = 1;
    while ~isempty(idx)
        idx = [];
        period_at_locs = 1./f_fit(locs/fs);
        diff_locs = diff(locs)/fs;
        for i=1:length(diff_locs)
            if diff_locs(i)<period_at_locs(i)/2
                idx = [idx;i];
            end
        end
        locs(idx+1)=[];
    end

    figure;
    plot(time, stim_smoothened); hold on;
    plot(locs/fs, stim_smoothened(locs), 'rx');
    yyaxis right; plot(cc/fs,rr,'k|');

    stim_clips = [];
    resp_clips = [];
    I_spike_phase = [];
    I_spike = [];
    spike_time_sd = [];
    fidelity = nan(1,length(locs));
    fid_freq = nan(1,length(locs));


    I = 1;
    for k= 3:length(locs)
        stim_clips = stim(1, locs(k-1):locs(k));
        resp_clips = resp(:, locs(k-1):locs(k));

        [row_ind,col_ind] = find(resp_clips == 1);


        [C, ia, ic] = unique(row_ind);
        ind = col_ind(ia);
        

        if length(col_ind)>=1
            fidelity(k) = length(ind)/P(end).complete_trials;
            fid_freq(k) = freq(ind(1)+locs(k-1));
            I_spike_time = ind/fs;
            if length(ind) >= 3
                spike_time_sd(I) = std(I_spike_time);
                I_spike(I) = freq(ind(1)+locs(k-1));%stim_freq(k-1);
                I=I+1;
            end
        end


    end

    figure;

    ax1 = subplot(3,1,1);
    sdfill(freq, P(end).avg_gcfr(startPt:stopPt), std(P(end).gcfr(:,startPt:stopPt), [],1),[0 0 0]);
    box off;
    axis padded;
    ax1.XAxis.Visible = "off";
    ylabel('Firing rate (Hz)');

    ax2 = subplot(3,1,2);
    scatter(I_spike, spike_time_sd*1000, '.', 'Color', 'k'); hold on;
    ylabel('Jitter (ms)');    
    yyaxis right; scatter(fid_freq, fidelity, '.', 'Color', 'k');
    ylabel('Fidelity');

    box off;
    axis padded;
    % ylim([0 8]);
    % yticks([0 4 8]);

    low_hz_ind = find(floor(I_spike)==60, 1,'first');
    mid_hz_low = find(floor(I_spike)==70,1,'first');
    mid_hz_high = find(floor(I_spike)==90,1,'first');
    high_hz_ind = find(floor(I_spike)==100,1,'first');

    jitter_1 = spike_time_sd(1:low_hz_ind);
    jitter_2 = spike_time_sd(mid_hz_low:mid_hz_high);
    jitter_3 = spike_time_sd(high_hz_ind:end);
    jitter_groups = [jitter_1'; jitter_2'; jitter_3'];

    g1 = repmat({'20-60 Hz'}, length(jitter_1),1);
    g2 = repmat({'70-90 Hz'}, length(jitter_2),1);
    g3 = repmat({'100-150 Hz'}, length(jitter_3),1);
    g=[g1;g2;g3];

    subplot(3,1,3);
    % boxplot(jitter_groups*1000, g)
    boxplot(jitter_groups*1000, g, 'Colors', 'k', 'Symbol','o', 'Jitter', 0.3);
    % ylabel('I spike jitter (ms)');
    % title('I spike jitter for frequency chirp stimulus')
    box off;
    axis padded;
    % ylim([0 8]);
    % yticks([0 4 8]);
    ylabel('Jitter (ms)');
    xlabel('Stimulus frequency (Hz)');


    h_j1 = kstest(jitter_1)
    h_j2 = kstest(jitter_2)
    h_j3 = kstest(jitter_3)

    ll = max([length(jitter_1) length(jitter_2) length(jitter_3)]);
    jitter_1(end+1:ll) = nan;
    jitter_2(end+1:ll) = nan;
    jitter_3(end+1:ll) = nan;
    [p,tbl,stats] = kruskalwallis([jitter_1' jitter_2' jitter_3'])

    multcompare(stats)
%}
    % pause;

end