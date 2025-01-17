function [pvalA, pvalAB, steadystateFR, baselineFR_mat, ss_position] = stepGCFRplots(P, c)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fs = P(1).fs;

labelFontSize = 12;
tickLabelSize = 10;
% c= parula(20);

stim_name = string(extractfield(P, 'stim_name'));
step_pos = str2double(extractAfter(stim_name, "amp_ "));
[~,idx] = sort(step_pos);
P = P(idx);


[b,a] = butter(3,4/(fs/2), 'low');
for i=1:length(P)
    P(i).mean_vel_filtered = filtfilt(b,a,(diff(P(i).mean_movement)).*fs);
    acceleration = diff(P(i).mean_vel_filtered)*fs;
    P(i).mean_acc_filtered = filtfilt(b, a, acceleration);
end

fig1 = plotPosVelAccGCFR(P, length(P), 0,0,0);
% savefigures(P(1), "traces", fig1, "png", 'D:\Work\Figures for presentation\step');
%Steady state firing rate Vs Position

% figure('WindowState','minimized');
% hold on;
steadystateFR = nan(P(1).no_of_trials, length(P));
position = [];
baselineFR_mat = nan(P(1).no_of_trials, length(P));
ss_position = nan(P(1).no_of_trials, length(P));
baseline_position = nan(P(1).no_of_trials, length(P));
% FRsorted = [];
% position_sorted = [];
meanFR = [];
% ssLoc = P(1).OFF_dur+3;

for i=1:length(P)
    stim_name = split(P(i).stim_name);
    amp = str2double(stim_name(2));

    % rampEndIdx = find(P(i).intendedStimulus(1,:)==amp,1, 'first');

    % ssLoc = rampEndIdx/P(1).fs + 1;
    [onLoc, offLoc] = miscFuncs.findSSbounds(P(i).mean_movement, 0.95, 10, P(i).fs);
    % ssBounds = [offLoc-2*P(1).fs  offLoc - 1*P(1).fs];
    ssBounds = [offLoc-1.5*P(1).fs  offLoc-0.5*P(i).fs];
    baselineBounds = [(P(i).OFF_dur-1.5)*fs+1 (P(i).OFF_dur-0.5)*fs];

    pos_ref = mean(P(i).antennal_movement(:,baselineBounds(1):baselineBounds(2)),2);
    pos_stim = mean(P(i).antennal_movement(:, ssBounds(1):ssBounds(2)),2);
    % pos_stim = (mean(P(i).antennal_movement(:, (ssLoc+1)*P(1).fs:(ssLoc+2)*P(1).fs),2));
    % pos_stim = (mean(P(i).antennal_movement(:, (ssLoc+2.5)*P(1).fs:(ssLoc+3.5)*P(1).fs),2));
    amplitude = pos_stim-pos_ref;
    position_values = repmat(mean(amplitude), [P(i).complete_trials,1]);
    baseline_pos = repmat(mean(pos_ref), [P(i).complete_trials,1]);

    ss_position(1:length(position_values),i) = position_values;
    baseline_position(1:length(baseline_pos),i) = baseline_pos;

    baselineFR =  mean(P(i).gcfr(:,baselineBounds(1):baselineBounds(2)),2);
    baselineFR_mat(1:length(baselineFR), i) = baselineFR;
    ssFR = mean(P(i).gcfr(:,ssBounds(1):ssBounds(2)),2);
    % ssFR = mean(P(i).gcfr(:,(ssLoc+1)*P(1).fs:(ssLoc+2)*P(1).fs),2); % earlier 1 s
    % ssFR = mean(P(i).gcfr(:,(ssLoc+2.5)*P(1).fs:(ssLoc+3.5)*P(1).fs),2);
    steadystateFR(1:length(ssFR), i) = ssFR;


end

pvalA = nan;
pvalAB = nan;


[xData, yData] = prepareCurveData([ss_position baseline_position], [steadystateFR baselineFR_mat]);
[fitresult, gof] = fit(xData, yData, 'poly1');

pBounds = predint(fitresult, xData, 0.95, 'functional', 'on');
% err = std(max_FR_sorted, [],1);
yfit = fitresult(xData);
rho = corr(xData, yData, 'type','Spearman')
if abs(rho) >=0.8
    figure;
    plot(xData,yData, 'Color', c,'Marker','.', 'MarkerSize',10, 'LineStyle','none'); hold on;
    plot(fitresult);
    % plot(fitresult, xData, yData, 'Color', c, 'LineStyle', '--');

    
    if fitresult.p1 >=0
    fill([unique(xData); flip(unique(xData))], [unique(pBounds(:,1)); flip(unique(pBounds(:,2)))],...
        c, 'FaceAlpha',0.2, 'EdgeColor','none');
    else
        fill([unique(xData); flip(unique(xData))], [flip(unique(pBounds(:,1))); (unique(pBounds(:,2)))],...
        c, 'FaceAlpha',0.2, 'EdgeColor','none');
    end
    

text(mean(xData),max(yData)-10, {string(gof.rsquare) string(fitresult.p1)});
% xticks([-1 -0.5 0 0.5 1]);
ylim([0 100]);
legend('', 'Box', 'off');
box off
end
% savefigures(P(1), 'linear_fit', gcf, "png", 'D:\Work\Figures for presentation\step\');

%{

meanFR = [mean(steadystateFR,1) mean(baselineFR_mat,"all")];
meanPos = [mean(ss_position,1) mean(baseline_position, "all")];








[posSorted, idx] = sort(meanPos, 'ascend');
FRsorted = meanFR(idx);


[pvalA, pvalAB] = miscFuncs.returnPvals(steadystateFR, baselineFR_mat);

% figure;
% if ~isempty(find(pvalA <= 0.01, 1)) % 2 stars
%     plot(posSorted, FRsorted, 'Color',c, 'LineStyle','-', 'Marker','.', 'MarkerSize',10); hold on;
% else
%     plot(posSorted, FRsorted, 'Color',[0 0 0 0.3], 'LineStyle','-', 'Marker','.', 'MarkerSize', 10); hold on;
% end
yshift = 1;
for i=2:length(pvalA)
    for j=1:i
        if pvalA(i,j)<0.01
            lineYLoc = max(steadystateFR(:,[i,j]),[],"all")+yshift;
            plot([ss_position(1,i),ss_position(1,j)],[lineYLoc lineYLoc], 'Color', [0 0 0 0.3]);
            text(mean([ss_position(1,i) ss_position(1,j)]),lineYLoc+1, "**", "FontSize", 20,"Color",[1 0 0 0.3]);
            yshift=yshift+1;
        end
    end
        end

for i=1:length(pvalAB)
    if pvalAB(i)<0.01
        lineYLoc = max([steadystateFR(:,i), baselineFR_mat(:,i)],[],"all")+yshift;
        plot([ss_position(1,i) baseline_position(1,i)],[lineYLoc lineYLoc], 'Color', [0 0 0 0.3]);
        text(mean([ss_position(1,i) baseline_position(1,i)]),lineYLoc+1, "**", "FontSize", 20,"Color",[1 0 0 0.3]);
        yshift = yshift+1;
    end
end

legend('', 'Box', 'off');
ylabel('Mean steady state firing rate (Hz)', 'FontSize',labelFontSize, 'Rotation',90);
xlabel('Position (deg)', 'FontSize',labelFontSize);
% title(replace([P(1).date P(1).filename], '_','-'));
ax = gca;
ax.FontSize = tickLabelSize;
ax.FontName = 'Calibri';
ax.Box = 'off';



% end
%}
end
