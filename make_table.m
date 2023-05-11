%% Table 
prompt = 'irow';
irow = input(prompt);
T.expt_data(irow) = expt_date;
T.filename(irow) = filename;
T.fs(irow) = fs;
T.STA(irow,:) = STA;
T.ev1(irow,:) = ev1';
T.ev2(irow,:) = ev2';
T.power_fft(irow,:) = power_fft;
T.frq_fft(irow,:) = frq_fft;
T.chirp_frq(irow,:) = inc_frq_chirp_f;
T.chirp_gcfr(irow,:) = mean(inc_chirp_gcfr,1);
T.velocity(irow,:) = nan;
T.maxFR(irow,:) = nan;
T.velocity(irow,1:length(vel)) = vel';
T.maxFR(irow,1:length(meanMaxFR)) = meanMaxFR;


%%
for irow=1:7
     [power_fft, frq_fft] = get_fft(T(irow,:).ev1, T.fs(irow), 0.03*T.fs(irow));  
     T.power_fft_ev1(irow,:) = power_fft;
     T.frq_fft_ev1(irow,:) = frq_fft;
end

%% Plots

figure;

for iplot = 1:6
    ax1 = subplot(2,6,iplot); plot(T.frq_fft(iplot,:), T.power_fft(iplot,:));
    title(replace(T.filename(iplot),'_',' '));
    xlim([0 150]);
    ax2 = subplot(2,6,iplot+6); plot(T.chirp_frq(iplot,:), T.chirp_gcfr(iplot,:));
    linkaxes([ax1 ax2], 'x');
end
