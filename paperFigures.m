%% Scripts used to generate paper figures

%% Fig 1B : Comparing generated and delivered stimulus

% 21.07.2022 M1 N3 used (this file does not exist)
% 12.07.2022 M1 N3

% miscFuncList = miscFuncs;
dataDirectory = '2022.06.12';
filename = "M1_N3"
gen_stim = stim.data.load;
gen_stim_1 = -gen_stim(:,1);


hes_stim_norm = miscFuncs.minmaxNorm(P(i).antennal_movement(1,:));
gen_stim_norm = miscFuncs.minmaxNorm(gen_stim_1');

figure;
plot(time(1:single_trial_length), gen_stim_norm); hold on;
plot(time(1:single_trial_length), hes_stim_norm);
xlim([2 14]);
box off;
ax = gca;
ax.YAxis.Visible = 'off';


%% Fig 1C : Basic data processing - 21.06.2022 M2 N5

% dataDirectory = '2022.07.07';
% filename = 'M1_N1_ramp'; % used in v1 of figures;

% Fig 1C v2 
% dataDirectory = '2022.07.13';
% filename = 'M1_N1_ramp';

% Fig 6A,B i (v3)
% 22.07.2022 M1_N1_stair
% 07.07.2022 M1_N1_stair

dataDirectory = "2022.07.22";
filename = "M1_N1_ramp";
P = getStructP(dataDirectory, filename,[nan nan],1);
%%
i=1;

figure;

tickLabelSize = 10;
labelFontSize = 12;

% stimulus plot
A4 = subplot(4,1,1);
patch([5,5,10,10],[-1,1,1,-1], [0.96,0.9,0.9]);
hold on;
sdfill(P(i).time(1:P(i).single_trial_length)', P(i).mean_movement, std(P(i).antennal_movement,[],1), [0.6, 0.2,0])
A4.Box = 'off';
A4.XAxis.Visible = 'off';


A4.FontSize = tickLabelSize;
A4.LineWidth =1;
ylabel('Position (deg)','FontSize', labelFontSize);
ylim([-1.5 1]);
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
ylabel('Membrane potential (mV)','FontSize', labelFontSize);
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
ylabel('Trials', 'FontSize', labelFontSize);
A2.Box = 'off';
A2.XAxis.Visible = 'off';
A2.FontName = 'Calibri';
%

% GCFR
A3 = subplot(4,1,4);
sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).gcfr,1), std(P(i).gcfr,[],1), [0.4660 0.6740 0.1880]); hold on;
yline(mean(P(i).avg_gcfr(1:4*P(i).fs)), 'k--');
ylim([-10 200])
A3.Box = 'off';
A3.XAxis.Visible = 'on';
A3.FontSize = tickLabelSize;
A3.LineWidth =1;
A3.FontName = 'Calibri';
ylabel('Firing rate (Hz)','FontSize', labelFontSize);
xlabel('Time(s)','FontSize', labelFontSize);


linkaxes([A1,A2,A3,A4], 'x');
xlim([4 13]);
% A3.XLim = [4 8];

%% Compare baseline and steady-state firing rate after stimulus

i=1;
[onLoc, offLoc] = miscFuncs.findSSbounds(P(i).mean_movement, 0.9, 10, P(i).fs);
ssBounds = [offLoc-1.5*P(i).fs  offLoc-0.5*P(i).fs];

figure; plot(1:length(P(i).mean_movement), [P(i).gcfr]); hold on; 
yyaxis right; plot(1:length(P(i).mean_movement), P(i).mean_movement);
xline([onLoc offLoc], 'k--'); 
xline(ssBounds, 'r--');
xline([(P(i).OFF_dur-1.5)*P(i).fs  (P(i).OFF_dur-0.5)*P(i).fs], 'r--');

% baselineMeanFR = mean(P(i).gcfr(:,3*fs:4*fs), 2);
% ssFR = mean(P(i).gcfr(:,6.5*fs:7.5*fs),2);

baselineMeanFR = mean(P(i).gcfr(:,(P(i).OFF_dur-1.5)*P(i).fs : (P(i).OFF_dur-0.5)*P(i).fs) , 2);
ssFR = mean(P(i).gcfr(:,ssBounds(1):ssBounds(2)),2);

figure;
scatter(ones([1,P(i).no_of_trials]), baselineMeanFR, 10,"black","filled");
hold on;
scatter(2*ones([1,P(i).no_of_trials]), ssFR, 10,"black","filled");
for i = 1:P(i).no_of_trials
    line([1,2], [baselineMeanFR(i),ssFR(i)]); hold on;
end
xlim([0 3]);
xticks(0:3);
xticklabels({'', 'Baseline position','1 deg constant deflection', ''})
ylabel('Firing rate (Hz)')
ylim([0 100])

pval = ranksum(baselineMeanFR, ssFR)
% pval = kstest2(baselineMeanFR, ssFR)

if pval <= 0.01
    text(1.5,max([baselineMeanFR; ssFR])-10, ["**"; pval])
end

%% Fig 2: Phasic and phaso-tonic responses
% 09.08.2022 M1_N2, M1_N4


dataDirectory = '2022.08.09';
filename = 'M1_N2_ramp';
P = getStructP(dataDirectory, filename,[nan nan],1);
%%
i=7;

[b,a] = butter(3,4/(P(i).fs/2), 'low');

tickLabelSize = 10;
labelFontSize = 12;

figure;
% Stimulus
A4 = subplot(4,1,1);
patch([P(i).OFF_dur,P(i).OFF_dur,P(i).OFF_dur+7,P(i).OFF_dur+7],[-1,0.2,0.2,-1], [0.96,0.9,0.9]); hold on;
sdfill(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, std(P(i).antennal_movement,[],1), [0.6, 0.2,0])

A4.Box = 'off';
A4.XAxis.Visible = 'off';
A4.FontSize = tickLabelSize;
A4.LineWidth =1;
ylabel('Position (deg)','FontSize', labelFontSize);
ylim([-1 0.2]);
% A4.FontName = 'Calibri';

% Velocity profile
A5 = subplot(4,1,2);
% velocity = diff(P(i).mean_movement)*fs;
velocity = diff(P(i).antennal_movement, 1, 2)*P(i).fs;
velocity = [velocity velocity(:,end)]; % Too make the length of velocity matrix same as time matrix
vel_filtered = filtfilt(b,a,velocity'); % Matrix transposed for filtfilt function
vel_filtered = vel_filtered';
% plot(time(2:single_trial_length),vel_filtered, 'LineWidth',2,'Color', [0.8824 0.6314 0], 'LineStyle',':');
sdfill(P(i).time(1:P(i).single_trial_length), mean(vel_filtered, 1), std(vel_filtered,[],1), [0.6, 0.2,0])
ylabel('Velocity (deg/s)', 'FontSize', labelFontSize);
ylim([-3.3 3])
yticks([-3 0 3])
A5.LineWidth =1;
A5.FontSize = tickLabelSize;
A5.Box = 'off';
A5.XAxis.Visible = 'off';
% A5.XLimitMethod = 'padded';
% A5.YLimitMethod = 'padded';
% A5.FontName = 'Calibri';

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
ylabel("Trials");
A2.XAxis.Visible = 'off';


% GCFR
A3 = subplot(4,1,4);
sdfill(P(i).time(1:P(i).single_trial_length), mean(P(i).gcfr,1), std(P(i).gcfr,[],1), [0.4660 0.6740 0.1880]); hold on;
yline(mean(P(i).avg_gcfr(1*P(i).fs : 3*P(i).fs)), 'k--');
[val,ind] = max(P(i).avg_gcfr);
% text(ind/P(i).fs, val+50, '\downarrow', 'Color','k' );
% text((ind/P(i).fs)+3, val+50, '\downarrow', 'Color','k' );
ylim([-10 200])
A3.Box = 'off';
% A3.XAxis.Visible = 'off';
A3.FontSize = tickLabelSize;
A3.LineWidth =1;
% A3.FontName = 'Calibri';
ylabel('Firing rate (Hz)','FontSize', labelFontSize);
xlabel('Time(s)','FontSize', labelFontSize);


linkaxes([A4, A5, A2, A3], 'x');
A4.XLim = [P(i).OFF_dur-2 P(i).OFF_dur+7];
%         linkaxes([A1,A2,A4], 'x');
% linkaxes([A3,A4,A5], 'x');

%% Range of adaptations

cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
rampData = readlines('ramp_forAdaptation.txt');
dataDirectory_ramp = extractBefore(rampData, '_');
filename_ramp = extractAfter(rampData, "_");
nfiles = length(filename_ramp);
gcfrSSdur = 3;
% ssStim = nan(1,nfiles);
% peakFR = nan(1,nfiles);
% baselineFR = nan(1,nfiles);
% ssFR = nan(1,nfiles);
varNames = {'dataDirectory' , 'filename', 'Stimulus', 'FiringRate', 'ssStim', 'peakFR', 'baselineFR', 'ssFR', 'Fs'};
T = table(varNames);

figure;
for irow = 1:nfiles
    expt_date = replace(dataDirectory_ramp(irow), '-','.');
    filenameParts = split(filename_ramp(irow), '_');
    mothId = filenameParts(1);
    if exist(join(['D:\Work\Recordings\' string(expt_date) '\' 'raw\' mothId '\'], ''), 'file') ~=0
        P = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),[nan nan],1);

        % for
        i = find(contains([P.stim_name], "1"));

        % stim_name = split(P(i).stim_name, " ");
        % delta_t = str2double(stim_name(2));
        delta_t = 0.25;

        if P(i).mean_movement((P(i).OFF_dur+1)*P(i).fs) < 0

        startPt = (P(i).OFF_dur - 1)*P(i).fs+1;
        stopPt = (P(i).OFF_dur + delta_t + gcfrSSdur)*P(i).fs;
        totalDur = 1 + delta_t + gcfrSSdur;

        % t = P(i).time(1:P(i).single_trial_length);
        % [~,maxFRLoc] = max(P(i).avg_gcfr);
        % stim = P(i).mean_movement(maxFRLoc:maxFRLoc+gcfr_dur*P(i).fs);
        % gcfr = P(i).avg_gcfr(maxFRLoc:maxFRLoc+gcfr_dur*P(i).fs);
        stim_norm = miscFuncs.minmaxNorm(P(i).mean_movement(startPt:stopPt));
        gcfr_norm = miscFuncs.minmaxNorm(P(i).avg_gcfr(startPt:stopPt));

        ssStim = P(i).mean_movement(stopPt);
        peakFR = max(P(i).avg_gcfr(startPt:stopPt));
        baselineFR = P(i).avg_gcfr(startPt);
        ssFR = P(i).avg_gcfr(stopPt);

        if P(i).fs ~=10000
            Stimulus = downsample(P(i).mean_movement(startPt:stopPt), P(i).fs/10000);
            FiringRate = downsample(P(i).avg_gcfr(startPt:stopPt), P(i).fs/10000);
        else
            Stimulus = P(i).mean_movement(startPt:stopPt);
            FiringRate = P(i).avg_gcfr(startPt:stopPt);
        end
        T(irow,varNames) = table(dataDirectory_ramp(irow), filename_ramp(irow), Stimulus, FiringRate, ssStim, peakFR, baselineFR, ssFR, P(i).fs);

        %{

            t = linspace(0, totalDur, length(Stimulus));
            ax1 = subplot(2,1,1);
            % plot(t, stim_norm, 'Color', [0.6, 0.2,0]); hold on;
            hold on;
            plot(t, stim_norm);
            % plot(t, P(i).mean_movement(startPt:stopPt)); %, 'Color', [0.6, 0.2,0]);
            ylabel('Angular position ({\circ})');
            ax1.XAxis.Visible = 'off';
            box off;

            ax2 = subplot(2,1,2);
            % plot(t, gcfr_norm, 'Color', [0.4660 0.6740 0.1880]); hold on;
            hold on;
            plot(t, gcfr_norm);
            
            % plot(t, P(i).avg_gcfr(startPt:stopPt)); %, 'Color',[0.4660 0.6740 0.1880]); hold on;
            % ylim([0 250]);
            ylabel('Firing rate (Hz)');
            xlabel('Time (s)');
            box off;
            linkaxes([ax1 ax2], 'x');
        %}
        % end
        end

    end
end

ind = find(T.ssStim==0);
T(ind,:)=[];

select_neurons = find(T.ssStim < median(T.ssStim)+0.05 & T.ssStim > median(T.ssStim)-0.05);

%%
idx = kmeans(T.Stimulus, 3);
c = [[0 0 0 0.8];...
     [0 0 0 0.5];...
     [1 0 0 0.8]];
totalDur = 6;
figure;
t = linspace(0, totalDur, length(T.FiringRate(1,:)));
for irow = 1:height(T) %select_neurons
    ax1 = subplot(3,1,1);
    hold on;
    plot(t, T.Stimulus(irow,:), 'Color', c(idx(irow),:)); %, 'Color', [0.6, 0.2,0])); hold on;
    % plot(t, P(i).mean_movement(startPt:stopPt)); %, 'Color', [0.6, 0.2,0]);
    ylabel('Angular position ({\circ})');
    ax1.XAxis.Visible = 'off';
    box off;
    
    ax2 = subplot(3,1,2);
    hold on;
    plot(t, T.FiringRate(irow,:), 'Color', c(idx(irow),:)); %, 'Color', [0.4660 0.6740 0.1880])); hold on;
    % plot(t, P(i).avg_gcfr(startPt:stopPt)); %, 'Color',[0.4660 0.6740 0.1880]); hold on;
    ylim([0 250]);
    % ylim([-0.1 1.1]);
    ylabel('Firing rate (Hz)');    
    box off;
    linkaxes([ax1 ax2], 'x');

    ax3 = subplot(3,1,3)
    hold on;
    plot(t, T.FiringRate(irow,:)-T.FiringRate(irow,1), 'Color', c(idx(irow),:)); %, 'Color', [0.4660 0.6740 0.1880])); hold on;
    % plot(t, P(i).avg_gcfr(startPt:stopPt)); %, 'Color',[0.4660 0.6740 0.1880]); hold on;
    ylim([0 250]);
    % ylim([-0.1 1.1]);
    ylabel({'Baseline subtracted'; 'Firing rate (Hz)'});
    xlabel('Time (s)');
    box off;
    linkaxes([ax1 ax2], 'x');

    
end


%% clustered data into subplots
ngrp = 3;
% idx = kmeans(T.FiringRate-T.FiringRate(:,1),ngrp);
[idx,C] = kmeans(T.ssStim, ngrp);
X=T.ssStim;
figure
gscatter(X(:,1),X(:,2),idx,'bgm')
hold on
plot(C(:,1),C(:,2),'kx')
legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid')

c = lines(ngrp);
% c = [[0 0 0 0.8]; [1 0 0 0.8]];
figure;
t = linspace(0, totalDur, length(T.FiringRate(1,:)));
for irow=1:height(T)
 
    ax1 = subplot(2,ngrp,idx(irow));
    hold on;
    plot(t, T.Stimulus(irow,:), 'Color', c(idx(irow),:)); %, 'Color', [0.6, 0.2,0])); hold on;
    % ylabel('Angular position ({\circ})');
    ax1.XAxis.Visible = 'off';
    ylim([-1.1 0.1]);
    box off;

    ax2 = subplot(2,ngrp,idx(irow)+ngrp);
    hold on;
    plot(t, T.FiringRate(irow,:)-T.FiringRate(irow,1), 'Color', c(idx(irow),:)); %, 'Color', [0.4660 0.6740 0.1880])); hold on;
    
    ylim([-50 200]);
    
    % ylabel({'Baseline subtracted'; 'Firing rate (Hz)'});
    xlabel('Time (s)');
    box off;
   
    linkaxes([ax1 ax2], 'x');

    
    % pause;

end
ssStim_grp4 = T.ssStim(idx==4)
std(ssStim_grp4)
ssStim_grp3 = T.ssStim(idx==3)
std(ssStim_grp3)
ssStim_grp2 = T.ssStim(idx==2)
std(ssStim_grp2)
ssStim_grp1 = T.ssStim(idx==1)
std(ssStim_grp1)


%% groups data by close values
totalDur = 6;
t = linspace(0, totalDur, length(T.FiringRate(1,:)));

figure; histogram(T.ssStim,'BinWidth', 0.02)
ylabel('Counts');
xlabel('Angular position ({\circ})');

% Choosing the middle group - between -0.82 and -0.72
idx = find(T.ssStim > -0.82 & T.ssStim < -0.72);
c = parula(height(T));
i=1;
for irow= idx
    c(irow,:)
    ax1 = subplot(2,1,1); 
    plot(t, T.Stimulus(idx,:), "Color" , [0.6, 0.2,0]); hold on;
    ylabel('Angular position ({\circ})');
    ax1.XAxis.Visible = 'off';
    ylim([-1.1 0.1]);
    box off;

    ax2 = subplot(2,1,2); hold on;
    plot(t, T.FiringRate(idx,:)-T.FiringRate(idx,1), 'Color', [0.4660 0.6740 0.1880]);
    ylim([-50 150]);    
    ylabel({'Baseline subtracted'; 'Firing rate (Hz)'});
    xlabel('Time (s)');
    box off;   
    linkaxes([ax1 ax2], 'x');
    i=i+1;
end

%% Quantifying adaptation
i=1;
t = P(i).time(1:P(i).single_trial_length);
[~,maxFRLoc] = max(P(i).avg_gcfr);
stim = P(i).mean_movement(maxFRLoc:maxFRLoc+3*P(i).fs);
gcfr = P(i).avg_gcfr(maxFRLoc:maxFRLoc+3*P(i).fs);
t1 = linspace((1/P(i).fs), 3, length(gcfr));
start_point = P(i).OFF_dur*P(i).fs+1;
stim_name_parts = split(P(i).stim_name); 
delta_t = str2double(stim_name_parts(2));
rampEndIdx = start_point+delta_t*P(i).fs - maxFRLoc;


[xData, yData] = prepareCurveData( t1', gcfr' );
% fo = fitoptions('Method', 'NonlinearLeastSquares', ...
%     'Robust', 'off');
% [fitresult, gof] = fit( xData, yData, 'power1', fo );
% tau1 = fitresult.b
% % tau2 = fitresult.d
% rsquare = gof.rsquare

model = @(b,xData) (b(1)*xData.^(-b(2)));
% model = @(b,xData) b(1)*exp(b(2)*xData);
beta0_empirical_phasic = [max(yData), 0.1];
beta0_refined_phasic = fminsearch(@(beta) norm(yData - model(beta, xData)), beta0_empirical_phasic);
mdl_phasic = fitnlm(xData, yData, model, beta0_refined_phasic);
y_pred = feval(mdl_phasic, t1); 
figure; plot(t1, gcfr, t1, y_pred);

% str = sprintf("{\\tau}_{1} = %0.3f \n{\\tau}_{2} = %0.3f \nrsquare = %0.3f", tau1, tau2, rsquare);
figure;
plot( fitresult, xData, yData ); hold on;
% ylim([0 200]);
% xlim([0 3]);
% yticks(0:40:200);
text(t1(end)-2,max(gcfr)-10, str);
ylabel('Mean firing rate (Hz)');
xlabel('Time (s)');
legend off;
box off;


% [xData, yData] = prepareCurveData( t1(1:rampEndIdx), gcfr(1:rampEndIdx) );
% fo = fitoptions('Method', 'NonlinearLeastSquares', ...
%     'Robust', 'off');
% [fitresult, gof] = fit( xData, yData, 'exp1', fo );
% tau1_separate = fitresult.b;
% plot(fitresult, 'c');
% xline(rampEndIdx/P(i).fs, 'k--');
% 
% [xData, yData] = prepareCurveData( t1(rampEndIdx+1:end), gcfr(rampEndIdx+1:end) );
% fo = fitoptions('Method', 'NonlinearLeastSquares', ...
%     'Robust', 'off');
% [fitresult, gof] = fit( xData, yData, 'exp1', fo );
% tau2_separate = fitresult.b;
% plot(fitresult, 'g');
% 



%% Fig 3 velocity encoding

% new - 14.07.2022 M1_N1_ramp

dataDirectory = '2022.07.14';
filename = 'M1_N1_ramp';
P = getStructP(dataDirectory, filename,[nan nan],1);
    
[velocities,amplitude_sorted, medMaxFR, k] = rampGCFRplots(P);
%%
figure;  hold on;
c = parula(height(T_vel_fr));
selectedCells = [];
for irow=1:height(T_vel_fr)
    [xData, yData] = prepareCurveData( cell2mat(T_vel_fr.velocity(irow,:)), cell2mat(T_vel_fr.fr(irow,:)));
    [fit_power, gof_power] = fit(xData, yData, 'power1');
    pBounds = predint(fit_power, xData, 0.95, 'functional', 'on');
    err = std(cell2mat(T_vel_fr.fr(irow)), [],1);
    yfit = fit_power(xData);
    % figure;
    if gof_power.rsquare >=0.8
        loglog(cell2mat(T_vel_fr.velocity(irow)), cell2mat(T_vel_fr.fr(irow)), 'Color', c(irow,:),'Marker','.', 'MarkerSize',10, 'LineStyle','none'); hold on;
        loglog(xData, yfit,'Color', c(irow,:));  
        fill([unique(xData); flip(unique(xData))], [unique(pBounds(:,1)); flip(unique(pBounds(:,2)))],...
            c(irow,:), 'FaceAlpha',0.2, 'EdgeColor','none');
        selectedCells = [selectedCells irow];
        % y_pred= feval(mdl, cell2mat(T_vel_fr.velocity(irow,:))); 
        % plot(cell2mat(T_vel_fr.velocity(irow,:)), y_pred); hold on;
        % errorbar(cell2mat(T_vel_fr.velocity(irow)), cell2mat(T_vel_fr.fr(irow)), err, 'ro',"MarkerSize",3,...
        % "MarkerEdgeColor","none","MarkerFaceColor",'k');
    end
    ax=gca;
    ax.XAxis.Scale = "log";
    ax.YAxis.Scale = "log";
    legend('', 'Box', 'off');
end
ylabel('Peak firing rate (Hz)');
xlabel('Angular velocity ({\circ}/s)');
xlim([0.1 5])
xticks([0.1 0.5 1 2 3 4 5])
ylim([25 250])
yticks([25 50 100 150 200 250])

%% Comparing steady-states for the same position following different velocities
stim_name = string(extractfield(P, 'stim_name'));
ramp_dur = str2double(extractAfter(stim_name, "ramp "));
[~,idx] = sort(ramp_dur);
P = P(idx);
n = numel(idx);

meanPos = [];
meanSSFR = [];
baselinePos = [];
baselineFR = [];
for i=1:n%length(P)
    [ss_onLoc, ss_offLoc] = miscFuncs.findSSbounds(abs(P(i).mean_movement), 0.95, 10, P(i).fs);
    % ssBounds = [ss_onLoc+0.5*P(i).fs ss_offLoc-0.5*P(i).fs];
    ssBounds = [ss_offLoc-1.5*P(i).fs  ss_offLoc-0.5*P(i).fs];
    baselineBounds = [(P(i).OFF_dur-1.5)*P(i).fs+1 (P(i).OFF_dur-0.5)*P(i).fs];
    
    meanPos = [meanPos mean(P(i).antennal_movement(:,ssBounds(1):ssBounds(2)), 2)];
    meanSSFR = [meanSSFR mean(P(i).gcfr(:,ssBounds(1):ssBounds(2)), 2)];
    baselinePos = [baselinePos mean(P(i).antennal_movement(:,baselineBounds(1):baselineBounds(2)), 2)];
    baselineFR = [baselineFR mean(P(i).gcfr(:, baselineBounds(1):baselineBounds(2)), 2)];

end

[pval,tbl,stats] = kruskalwallis(meanSSFR);
multcompare(stats)
[pval_bss,tbl_bss,stats_bss] = kruskalwallis([baselineFR meanSSFR]);
multcompare(stats_bss)

c = [0.1765    0.1804    0.5137;
    0.1765    0.1804    0.5137;
    0.1765    0.1804    0.5137;
    0.9137    0.3059    0.1059;
    0.9137    0.3059    0.1059];
alphaVal  = [1 0.7 0.5 0.7 1];
figure;
for i=1:n
    scatter([baselinePos(:,i) meanPos(:,i)], [baselineFR(:,i) meanSSFR(:,i)], [], "filled", 'MarkerFaceColor', c(i,:), "MarkerFaceAlpha",alphaVal(i)); hold on;
end

xticks([-0.8 0]);
yticks(50:10:100)
ylim([50 100])
ax=gca;
ax.XAxis.Direction = "reverse";
ylabel("Steady state firing rate (Hz)");
xlabel("Angular position ({\circ})");

%% Fig 4 Position encoding neurons
% dataDirectory = '2022.07.18';
% filename = 'M1_N4_step';

close all
dataDirectory = '2022.08.24';
filename = 'M1_N4_ramp';
filename_ramp = 'M1_N4_ramp';

rampGCFRplots(P)

%% Unusual neurons
% 2022.07.12_M1_N3_step
% 2022.07.12_M1_N3_ramp
% 18.07.2022 M1_N4_step and _ramp

dataDirectory = '2022.08.24';
filename = 'M1_N4_step';
P = getStructP(dataDirectory, filename,[nan nan],1);
stim_name = string(extractfield(P, 'stim_name'));
ramp_dur = str2double(extractAfter(stim_name, "amp_ "));
[~,idx] = sort(ramp_dur);
P = P(idx);
% n = 4%numel(idx);
% rampGCFRplots(P)
% protocolPlot(P);
[pvalA, pvalAB, steadystateFR, baselineFR_mat, ss_pos] = stepGCFRplots(P);

max_gcfr_sub = [];
% Peak firing rate in step protocol
P_sub = P;
ss_pos_sub = ss_pos;
% P_sub = P([4 5 6]);
% ss_pos_sub  = ss_pos(:,[4 5 6]);
for i = 1:length(P_sub)
    max_gcfr_sub = [max_gcfr_sub max(P_sub(i).gcfr,[],2)];
end

figure;
scatter(ss_pos_sub, max_gcfr_sub);

[p,tbl,stats] = kruskalwallis(max_gcfr_sub)
multcompare(stats)

%% Fig 5 Stair protocol
% 01.07.2022 M1_N1_stair2 %old

% 22.07.2022 M1_N1_stair
% 07.07.2022 M1_N1_stair

dataDirectory_stair = '2022.07.07';
filename_stair = 'M1_N1_stair'; 
downsampleFactor = 2;
P_stair = getStructP(dataDirectory_stair, filename_stair,[nan nan], downsampleFactor);

stairGCFRplots(P_stair);
% exportgraphics(figHandle, 'Fig5Ai-bottom.eps', 'ContentType', 'vector');


%% Adaptation figure
% 18.07 M1N4
% 07.07 M2N3
% 07.07 M1N1

dataDirectory = '2022.07.12';
filename = 'M1_N3_step';
P = getStructP(dataDirectory, filename,[nan nan],1);
% stim_name = string(extractfield(P, 'stim_name'));
% ramp_dur = str2double(extractAfter(stim_name, "ramp "));
% [~,idx] = sort(ramp_dur);
% P = P(idx);
% n = 4%numel(idx);
newColors = [0 0 0 1;
            0 0 0 0.5;
            0.8 0 0 0.5;
            0.8 0 0 1];
% rampGCFRplots(P)
% figure; protocolPlot(P); hold on;
for i=1:4%length(P)
    % [~, ~,~, ~, ax] = estimatePhasicAdaptation(P(i), newColors(i,:)); hold on;

    estimateTonicAdaptation(P(i), newColors(i,:));
end
% hold on;
% gca = ax;
% estimateTonicAdaptation(P(1));

% uniqueGroups = unique(T_phasic.neuronID);
% 

