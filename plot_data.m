function plot_data(P, varargin)
% Valid varargin: stimulus, velocity, membrane potential, raster, gcfr
waveforms = string(varargin)
disp(nargin)

nplots = nargin-1
isubplot = 0;

for i=1:length(P)

    if isempty(P(i).gcfr)
        continue;
    end

    fig = figure;

    if sum(contains(waveforms, "stimulus"))
        isubplot = isubplot+1;
        a(isubplot) = subplot(nplots,1,isubplot);
        % sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).stim_ifb,1), std(P(i).stim_ifb,[],1), [0.6, 0.2,0], "none")
        sdfill(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, std(P(i).antennal_movement,[],1), [0.6, 0.2,0], "")
        a(isubplot).Box = 'off';
        a(isubplot).XAxis.Visible = 'off';
        ylabel('Position (deg)');
        % ylim([-1 1]);
    end

    if sum(contains(waveforms, "velocity"))
        [x,y] = butter(3,4/(P(i).fs/2), 'low');
        isubplot = isubplot+1;
        a(isubplot) = subplot(nplots,1,isubplot);
        velocity = diff(P(i).mean_movement)*P(i).fs;
        vel_filtered = filtfilt(x,y,velocity);
        plot(P(i).time(2:P(i).single_trial_length),vel_filtered, 'LineWidth',2,'Color', [0.8824 0.6314 0], 'LineStyle',':');
        ylabel('Velocity (deg/s)');
        ylim([-2 2])
        a(isubplot).Box = 'off';
        a(isubplot).XAxis.Visible = 'off';
    end

    if sum(contains(waveforms, "membrane potential"))
        isubplot = isubplot+1;
        [p,l] = findpeaks(P(i).rec(1,:), "MinPeakHeight",0.25*max(P(i).rec(1,:)));
        a(isubplot) = subplot(nplots,1,isubplot);
        plot(P(i).time(1:P(i).single_trial_length), P(i).rec(1, :),'LineWidth', 0.5, 'Color', 'k'); %#0072BD');
        hold on; plot((l/P(i).fs), p, '.', 'MarkerEdgeColor', 'r', 'MarkerSize', 8); %'#A2142F');
        a(isubplot).Box = 'off';
        a(isubplot).XAxis.Visible = 'off';
        ylabel('Membrane potential (mV)');
    end

    if sum(contains(waveforms, "raster"))
        isubplot = isubplot+1;
        a(isubplot) = subplot(nplots,1,isubplot);

        k = 0.5;
        for j = 1:P(i).complete_trials
            l = find(P(i).raster(j,:)==1);
            spike_P(i).time = l/P(i).fs;
            for m = 1:length(spike_P(i).time)
                line([spike_P(i).time(m) spike_P(i).time(m)], [k k+0.5], 'Color', 'k', 'LineWidth', 1);
            end
            k = k+1;
        end

        ylabel('Trials');%, 'FontSize', labelFontSize, 'Rotation', 0);
        a(isubplot).Box = 'off';
        a(isubplot).XAxis.Visible = 'off';
    end


    if sum(contains(waveforms, "gcfr"))

        isubplot = isubplot + 1;
        a(isubplot) = subplot(nplots,1,isubplot);
        sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).gcfr,1), std(P(i).gcfr,[],1), [0.4660 0.6740 0.1880], "none") 
        ylabel('Firing rate (Hz)');
        xlabel('Time (s)');
        % ylim([0 200])
        a(isubplot).Box = 'off';
        a(isubplot).XAxis.Visible = 'off';
    end

    linkaxes(a, 'x');

    % if (P(i).stim_name == "sin" || P(i).stim_name == "sqr" || P(i).stim_name == "sum_sine")
    %     title_str = ((join(split(P(i).stim_name,"_")," ")) +" Hz");
    % else
    %     title_str = ((join(split(P(i).stim_name,"_"))));
    % end

    title_str = replace([string(P(i).date) P(i).filename], "_", " ");

    a(1).Title.String = title_str;
    a(end).XAxis.Visible = 'on';

end

end