%% Scripts used to generate paper figures

%% Fig 1B : Comparing generated and delivered stimulus

% 21.07.2022 M1 N3 used (this file does not exist)
miscFuncList = miscFuncs;

gen_stim = stim.data.load;
gen_stim_1 = -gen_stim(:,1);


hes_stim_norm = miscFuncList.minmaxNorm(P(1).antennal_movement(1,:));
gen_stim_norm = miscFuncList.minmaxNorm(gen_stim_1');

figure;
plot(time(1:single_trial_length), gen_stim_norm); hold on;
plot(time(1:single_trial_length), hes_stim_norm);
xlim([2 14]);
box off;
ax = gca;
ax.YAxis.Visible = 'off';


%%
%Fig 1C : Basic data processing - 21.06.2022 M2 N5
i=1 ;


tickLabelSize = 10;
labelFontSize = 12;

% stimulus plot
A4 = subplot(4,1,1); 
patch([5,5,8,8],[-1,1,1,-1], [0.96,0.9,0.9]);
hold on;
sdfill(P(i).time(1:P(i).single_trial_length), -P(i).mean_movement, std(-P(i).antennal_movement,[],1), [0.6, 0.2,0])
A4.Box = 'off';
A4.XAxis.Visible = 'off';


A4.FontSize = tickLabelSize;
A4.LineWidth =1;
ylabel('Position (deg)','FontSize', labelFontSize, 'Rotation', 0);
ylim([-1 1]);
A4.FontName = 'Calibri';

% Single trial membrane potential
[p,l] = findpeaks(P(i).rec(1,:), "MinPeakHeight",0.25*max(P(i).rec(1,:)));
A1 = subplot(4,1,2);
plot(P(i).time(1:P(i).single_trial_length), P(i).rec(1, :),'LineWidth', 0.5, 'Color', 'k'); %#0072BD');
hold on; plot((l/P(i).fs), p, '.', 'MarkerEdgeColor', 'r', 'MarkerSize', 8); %'#A2142F');
hold off;
ylim([-20 60]);
A1.Box = 'off';
A1.XAxis.Visible = 'off';
A1.FontSize = tickLabelSize;
A1.LineWidth = 1;
ylabel('Membrane potential (mV)','FontSize', labelFontSize, 'Rotation', 0);
A1.FontName = 'Calibri';


% raster plot
A2 = subplot(4,1,3);
k = 0.5;
for j = 1:P(i).complete_trials
    l = find(P(i).raster(j,:)==1);
    spike_time = l/P(i).fs;
    for m = 1:length(spike_time)
        line([spike_time(m) spike_time(m)], [k k+0.5], 'Color', 'k', 'LineWidth', 1);
    end
    k = k+1;
end
A2.FontSize = tickLabelSize;
A2.LineWidth =1;
A2.YLimitMethod = 'padded';
ylabel('Trials', 'FontSize', labelFontSize, 'Rotation', 0);
A2.Box = 'off';
A2.XAxis.Visible = 'off';
A2.FontName = 'Calibri';
%

% GCFR
A3 = subplot(4,1,4); 
sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).gcfr,1), std(P(i).gcfr,[],1), [0.4660 0.6740 0.1880]); hold on;
yline(mean(P(i).avg_gcfr(1:4*fs)), 'k--');
ylim([0 150])
A3.Box = 'off';
A3.XAxis.Visible = 'on';
A3.FontSize = tickLabelSize;
A3.LineWidth =1;
A3.FontName = 'Calibri';
ylabel('Firing rate (Hz)','FontSize', labelFontSize, 'Rotation', 0);
xlabel('Time(s)','FontSize', labelFontSize);


linkaxes([A1,A2,A3,A4], 'x');
A3.XLim = [4 8];

%% Compare baseline and steady-state firing rate after stimulus

baselineMeanFR = mean(P(1).gcfr(:,3*fs:4*fs), 2);
ssFR = mean(P(1).gcfr(:,6.5*fs:7.5*fs),2);

figure;
scatter(ones([1,P(1).no_of_trials]), baselineMeanFR, 10,"black","filled");
hold on;
scatter(2*ones([1,P(1).no_of_trials]), ssFR, 10,"black","filled");
for i = 1:P(1).no_of_trials
    line([1,2], [baselineMeanFR(i),ssFR(i)]); hold on;
end
xlim([0 3]);
xticks(0:3);
xticklabels({'', 'Baseline position','1 deg constant deflection', ''})
ylabel('Firing rate (Hz)')
ylim([0 50])

pval = ranksum(baselineMeanFR, ssFR)

if pval < 0.01
    text(1.5,max([baselineMeanFR; ssFR])-10, "*")
end

%% Fig 2: Phasic and phaso-tonic responses
% 09.08.2022 M1_N2, M1_N4



i=6 ;

[b,a] = butter(3,4/(fs/2), 'low');

tickLabelSize = 10;
labelFontSize = 12;

figure;
% Stimulus
A4 = subplot(3,1,1); 
patch([5,5,12,12],[-1,0.2,0.2,-1], [0.96,0.9,0.9]);
sdfill(P(i).time(1:P(i).single_trial_length), -P(i).mean_movement, std(-P(i).antennal_movement,[],1), [0.6, 0.2,0])

A4.Box = 'off';
A4.XAxis.Visible = 'off';
A4.FontSize = tickLabelSize;
A4.LineWidth =1;
ylabel('Position (deg)','FontSize', labelFontSize, 'Rotation', 0);
ylim([-1 0.2]);
A4.FontName = 'Calibri';

% Velocity profile
A5 = subplot(3,1,2);
% velocity = diff(-P(i).mean_movement)*fs;
velocity = diff(P(i).antennal_movement, 1, 2)*fs;
velocity = [velocity velocity(:,end)];
vel_filtered = filtfilt(b,a,velocity');
% plot(time(2:single_trial_length),vel_filtered, 'LineWidth',2,'Color', [0.8824 0.6314 0], 'LineStyle',':');
sdfill(P(i).time(1:P(i).single_trial_length), mean(vel_filtered, 2), std(vel_filtered,[],2), [0.6, 0.2,0])
ylabel('Velocity (deg/s)', 'FontSize', labelFontSize, 'Rotation', 0);
ylim([-3 3])
yticks([-3 0 3])
A5.LineWidth =1;
A5.FontSize = tickLabelSize;
A5.Box = 'off';
A5.XAxis.Visible = 'off';
% A5.XLimitMethod = 'padded';
% A5.YLimitMethod = 'padded';
A5.FontName = 'Calibri';

% GCFR
A3 = subplot(3,1,3); 
sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).gcfr,1), std(P(i).gcfr,[],1), [0.4660 0.6740 0.1880])
ylim([-10 150])
A3.Box = 'off';
A3.XAxis.Visible = 'on';
A3.FontSize = tickLabelSize;
A3.LineWidth =1;
A3.FontName = 'Calibri';
ylabel('Firing rate (Hz)','FontSize', labelFontSize, 'Rotation', 0);
xlabel('Time(s)','FontSize', labelFontSize);


linkaxes([A4, A5, A3], 'x');
A3.XLim = [3 12];
%         linkaxes([A1,A2,A4], 'x');
% linkaxes([A3,A4,A5], 'x');


%% Fig 3 velocity and position encoding & Fig 4
% 13.07.20222 M1_N1_ramp
% new - 14.07.2022 M1_N1_ramp
% 18.07.2022 M1_N4



% rampGCFRplots(P,c)
[pvalA, pvalAB, steadystateFR, baselineFR_mat] = stepGCFRplots(P);


%% Fig 5 Stair protocol
% 01.07.2022 M1_N1_stair2
% 07.07.2022 M1_N1_stair

dataDirectory_stair = '2022.07.01';
filename_stair = 'M1_N1_stair2';
downsampleFactor = 2;
P_stair = getStructP(dataDirectory_stair, filename_stair,0, downsampleFactor);
%%

figHandle = stairGCFRplots(P_stair);
exportgraphics(figHandle, 'Fig5Ai-bottom.eps', 'ContentType', 'vector');


