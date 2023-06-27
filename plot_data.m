function plot_data(single_trial_length,no_of_protocols, fs, time, filename,  P)

    [b,a] = butter(3,4/(fs/2), 'low');
    labelFontSize = 14;
    tickLabelSize = 12;
    
    
    for i=1:no_of_protocols
        
        if isempty(P(i).gcfr)
            continue;
        end
        
        fig = figure;
        

        % [p,l] = findpeaks(P(i).rec(1,:), "MinPeakHeight",0.25*max(P(i).rec(1,:)));
        % A1 = subplot(4,1,2);
        % plot(time(1:single_trial_length), P(i).rec(1, :),'LineWidth', 0.5, 'Color', 'k'); %#0072BD');
        % hold on; plot((l/fs), p, '.', 'MarkerEdgeColor', 'r', 'MarkerSize', 8); %'#A2142F');
        % hold off;
        % A1.Box = 'off';
        % A1.XAxis.Visible = 'off';
        % A1.FontSize = tickLabelSize;
        % A1.LineWidth = 1;
        % ylabel('Membrane potential (mV)','FontSize', labelFontSize, 'Rotation', 0);
        % A1.FontName = 'Calibri';
% 
%         A2 = subplot(4,1,3);
% % %         A2 = subplot(2,1,2);
% % 
%         k = 0.5;
%         for j = 1:P(i).complete_trials
%             l = find(P(i).raster(j,:)==1);
%             spike_time = l/fs;
%             for m = 1:length(spike_time)
%                 line([spike_time(m) spike_time(m)], [k k+0.5], 'Color', 'k', 'LineWidth', 1);
%             end
%             k = k+1;
%         end
%         A2.FontSize = tickLabelSize;
%         A2.LineWidth =1;
%         A2.XLimitMethod = 'padded';
%         A2.YLimitMethod = 'padded';
%         ylabel('Trials', 'FontSize', labelFontSize, 'Rotation', 0);
%         A2.Box = 'off';
%         A2.XAxis.Visible = 'off';
        % A2.FontName = 'Calibri';
% 
        
            
        % A3 = subplot(4,1,4); %plot(time(1:single_trial_length), P(i).norm_gcfr, 'Color', [0.2,0.3,0.49]);
        A3 = subplot(3,1,3);
%           % A3 = subplot(2,1,2);
        % A3 = subplot(2,1,2); 
        % [lineOut, ~] = stdshade(P(i).gcfr,0.2,[0.4660 0.6740 0.1880],time(1:single_trial_length)); 
        % lineOut.LineWidth  = 1;
        sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).gcfr,1), std(P(i).gcfr,[],1), [0.4660 0.6740 0.1880])
        ylabel('Firing rate (Hz)','FontSize', labelFontSize, 'Rotation', 0);
        xlabel('Time(s)','FontSize', labelFontSize);
        ylim([0 200])
        A3.Box = 'off';
        A3.XAxis.Visible = 'on';
        A3.FontSize = tickLabelSize;
        A3.LineWidth =1;
        % A3.XLimitMethod = 'padded';
        % A3.YLimitMethod = 'padded';
        A3.FontName = 'Calibri';
        
        % A4 = subplot(2,1,1);
        A4 = subplot(3,1,1);
        % A4 = subplot(4,1,1); %plot(time(1:single_trial_length), mean(P(i).antennal_movement), 'Color', [0.6, 0.2,0]);
        % [lineOut, ~] = stdshade(-P(i).antennal_movement,0.2,[0.6, 0.2,0],time(1:single_trial_length));
        % lineOut.LineWidth = 1;
        sdfill(P(i).time(1:P(i).single_trial_length), -P(i).mean_movement, std(-P(i).antennal_movement,[],1), [0.6, 0.2,0])
        A4.Box = 'off';
        A4.XAxis.Visible = 'off';
        A4.FontSize = tickLabelSize;
        A4.LineWidth =1;
%         ylabel('Indenter feedback voltage');                     
        ylabel('Position (deg)','FontSize', labelFontSize, 'Rotation', 0);
        ylim([-1 1]);
        % A4.XLimitMethod = 'padded';
        % A4.YLimitMethod = 'padded';
        A4.FontName = 'Calibri';

        A5 = subplot(3,1,2);
        velocity = diff(-P(i).mean_movement)*fs;
        vel_filtered = filtfilt(b,a,velocity);
        plot(time(2:single_trial_length),vel_filtered, 'LineWidth',2,'Color', [0.8824 0.6314 0], 'LineStyle',':');
        ylabel('Velocity (deg/s)', 'FontSize', labelFontSize, 'Rotation', 0);
        ylim([-2 2])
        A5.LineWidth =1;
        A5.FontSize = tickLabelSize;
        A5.Box = 'off';
        A5.XAxis.Visible = 'off';
        % A5.XLimitMethod = 'padded';
        % A5.YLimitMethod = 'padded';
        A5.FontName = 'Calibri'; 
        % if (P(i).stim_name == "sin" || P(i).stim_name == "sqr" || P(i).stim_name == "sum_sine")
        %     title((join(split(P(i).stim_name,"_")," ")) +" Hz");
        % else
        %     title((join(split(P(i).stim_name,"_"))));
        % end
        
        % linkaxes([A1,A2,A3,A4], 'x');
%         linkaxes([A1,A2,A4], 'x');
          linkaxes([A3,A4,A5], 'x');
    
%         savefigures(filename, P(i).stim_name, fig, P(i).date);
    end
    
end    