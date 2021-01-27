for i=1:9
t = linspace(0,16, length(T_N6.avg_gcfr(1,:)));
figure;
    for j=1:5
        a1 = subplot(2,1,1); plot(t, T_N6.gcfr{i}(j,:)/max(T_N6.gcfr{i}(j,:))); hold on;
        
        ylabel 'Avg GCFR';
        title (strcat("stim period", "=", num2str(T_N6.stim_period(i)), " s"));
        a2 = subplot(2,1,2); plot(t, T_N6.antennal_movement{i}(j,:)); hold on;
        ylabel 'Antennal movement'
        xlabel 'time (s)';
        linkaxes([a1 a2], 'x');
    end

end