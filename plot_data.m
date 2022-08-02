function plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P)

    for i=1:no_of_protocols
        
        if isempty(P(i).gcfr)
            continue;
        end
        
        fig = figure;
%         fig = figure(i+1);
        fig.Position = [0 0 960 640];
%         [p,l] = findpeaks(P(i).rec(1,:), "MinPeakHeight",0.25*max(P(i).rec(1,:)));
%         A1 = subplot(4,1,2);
%         plot(time(1:single_trial_length), P(i).rec(1, :),'LineWidth', 0.01, 'Color', 'k'); %#0072BD');
%         hold on; plot((l/fs), p, '.', 'MarkerEdgeColor', 'r'); %'#A2142F');
%         hold off;
%         A1.Box = 'off';
%         A1.XAxis.Visible = 'off';
%         ylabel('Membrane potential (mV)');
%         
        
%         A2 = subplot(4,1,3);
        A2 = subplot(3,1,2);
        
        k = 0.5;
        for j = 1:P(i).complete_trials
            l = find(P(i).raster(j,:)==1);
            spike_time = l/fs;
            for m = 1:length(spike_time)
                line([spike_time(m) spike_time(m)], [k k+0.5], 'Color', 'k', 'LineWidth', 0.5);
            end
            k = k+1;
        end
        ylabel('Trials', 'FontSize', 14);
        A2.Box = 'off';
        A2.XAxis.Visible = 'off';
        A2.YAxis.FontSize = 12;
        
            
%         A3 = subplot(4,1,3); plot(time(1:single_trial_length), P(i).norm_gcfr, 'Color', [0.2,0.3,0.49]);
        A3 = subplot(3,1,3);
%           A3 = subplot(2,1,2);
%         A3 = subplot(4,1,4); 
        [lineOut, ~] = stdshade(P(i).gcfr,0.2,[0.4660 0.6740 0.1880],time(1:single_trial_length)); %10 = (fs/L)*gcfr Hz
        lineOut.LineWidth  = 0.05;
        ylabel('Firing rate (Hz)','FontSize', 14);
        xlabel('time(s)','FontSize', 14);
        A3.Box = 'off';
        A3.XAxis.Visible = 'on';
        A3.YAxis.FontSize = 12;
        A3.XAxis.FontSize = 12;
        
%         A4 = subplot(2,1,1);
        A4 = subplot(3,1,1);
%         A4 = subplot(4,1,1); %plot(time(1:single_trial_length), mean(P(i).antennal_movement), 'Color', [0.6, 0.2,0]);
        [lineOut, ~] = stdshade(-P(i).antennal_movement,0.2,[0.6, 0.2,0],time(1:single_trial_length));
        
        lineOut.LineWidth = 0.05;
        A4.Box = 'off';
        A4.XAxis.Visible = 'off';
%         ylabel('Indenter feedback voltage');
        ylabel('Stimulus (mm)','FontSize', 14);
        A4.YAxis.FontSize = 12;
        
        if (P(i).stim_name == "sin" || P(i).stim_name == "sqr" || P(i).stim_name == "sum_sine")
            title((join(split(P(i).stim_name,"_")," ")) +" Hz");
        else
            title((join(split(P(i).stim_name,"_"))));
        end
        
%         linkaxes([A1,A2,A3,A4], 'x');
        linkaxes([A2,A3,A4], 'x');
%           linkaxes([A3,A4], 'x');
    
%         savefigures(filename, P(i).stim_name, fig, P(i).date);
    end
    
end    