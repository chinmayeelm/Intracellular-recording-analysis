for i=1:11
    % diff_sd = max_sd - min_sd;
    figure(1);
    hold on;
    plot( S(i).max_stim_sd, S(i).jitter, '.'); 
    xlabel('Max STD in stimulus patterns (m)');
    ylabel('Jitter (s)'); 
    % title('58');
    
    
    figure(2);
    hold on;
%     colormap turbo(11);
    histogram(S(i).jitter,'BinWidth', 1000/S(i).fs); 
    xlabel('Jitter (ms)');
    ylabel('Occurances');
    xlim ([0 2]);
    
    % figure;
    % plot(freq, S(i).jitter, '.');
    % xlabel('Max frequency in spike triggering stimulus pattern (Hz)');
    % ylabel('Jitter (ms)');
    
    figure(3);
    hold on;
    plot(S(i).slope, S(i).jitter,  '.'); 
    xlabel('Max amplitude in spike triggering stimulus pattern (mm)');
    ylabel('Jitter (ms)'); 
    
    figure(4);
    hold on;
    plot(S(i).pos_amp, S(i).jitter, '.'); 
    xlabel('Slope of spike triggering stimulus pattern (mm/sec)');
    ylabel('Jitter (ms)');
    
    figure(5);
    hold on;
    plot(S(i).stim_pattern_r, S(i).jitter, '.'); 
    xlabel('Correlation with STA');
    ylabel('Jitter (ms)');
    
    figure(6);
    hold on;
    plot(S(i).time_to_prev_spike, S(i).jitter, '.'); 
    xlabel('Time to previous spike (ms)');
    ylabel('Jitter (ms)'); 

%     colormap turbo(11);
    figure(7);
    hold on;
    scatter(S(i).baseline_FR, mean(S(i).jitter), 40, 'o', 'filled'); 
    xlabel('Baseline firing rate (Hz)');
    ylabel('Jitter (ms)'); 
    ylim ([0 2]);
    
%     colormap turbo(11);
    figure(8);
    hold on;
    scatter(diff(S(i).time_to_prev_spike_column), S(i).jitter(2:end), 100, '.'); 
    xlabel('Time to previous spiking event (ms)');
    ylabel('Jitter (ms)'); 
    ylim ([0 2]);
    
end
