function P = plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P)

    for i=1:no_of_protocols
        
        fig = figure(i);
        [p,l] = findpeaks(P(i).rec(1,:), "MinPeakHeight",0.2*max(P(i).rec(1,:)));
        A1 = subplot(4,1,1); plot(time(1:single_trial_length), P(i).rec(1, :),'LineWidth', 0.01, 'Color', '#0072BD');
        hold on; plot(l/fs, p, '.', 'MarkerEdgeColor', '#A2142F');
        ylabel('Voltage (mV)');
        title((join(split(P(i).stim_name,"_")," ")) +" Hz");
        
        A2 = subplot(4,1,2);
        k = 0.5;
        for j = 1:P(i).complete_trials
            l = find(P(i).raster(j,:)==1);
            spike_time = l/fs;
            for m = 1:length(spike_time)
                line([spike_time(m) spike_time(m)], [k k+0.5], 'Color', 'k', 'LineWidth', 0.01);
            end
            k = k+1;
        end
        ylabel('Trials');
            
%         A3 = subplot(4,1,3); plot(time(1:single_trial_length), P(i).norm_gcfr, 'Color', [0.2,0.3,0.49]);
        A3 = subplot(4,1,3); [lineOut, ~] = stdshade(P(i).norm_gcfr,0.6,[0.4660 0.6740 0.1880],time(1:single_trial_length));
        lineOut.LineWidth  = 0.05;
        ylabel('Normalised GCFR');
        
        
        A4 = subplot(4,1,4); %plot(time(1:single_trial_length), mean(P(i).antennal_movement), 'Color', [0.6, 0.2,0]);
        [lineOut, ~] = stdshade(P(i).antennal_movement,0.6,[0.6, 0.2,0],time(1:single_trial_length));
        lineOut.LineWidth = 0.05;
        
%         ylabel('Indenter feedback voltage');
        ylabel('Antennal movement');
        xlabel('time(s)');
        
        linkaxes([A1,A2,A3,A4], 'x');
    
%         savefigures(filename, P(i).stim_name, fig, P(i).date);
    end
    
end    