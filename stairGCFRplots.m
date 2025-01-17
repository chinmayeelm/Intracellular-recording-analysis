function hysteresis_area = stairGCFRplots(P)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
total_trial_dur = P.single_trial_length/P.fs;
time = linspace(0,total_trial_dur,P(1).single_trial_length);

stimulus = P.mean_movement;
mid_point = length(stimulus)/2;
labelFontSize = 11;
tickLabelSize = 10;

S = abs(diff(P.intendedStimulus(1,1:P.single_trial_length)));

x = [];
for i=1:length(S)-1
    if (S(i)==0 && S(i+1)>0) || (S(i)>0 && S(i+1)==0)
        x = [x i];
    end
end
x = reshape(x,2,[]);
x(3,:) = x(2,:);
x(4,:) = x(1,:);
x = x/P.fs;

% figure;
% plot(P.time(1:P.single_trial_length), P.mean_movement); hold on;
% xline(x(1,:), 'k--');

y_stim = repmat([-1.2 -1.2 1.2 1.2]', 1,length(x));
y_gcfr = repmat([-10 -10 150 150]', 1,length(x));


% Stimulus and response
%{
figHandle = figure('Color', 'w', 'Visible','on');%, 'WindowState','minimized');
ax1 = subplot(2,1,1);
plot(P.time(1:P.single_trial_length), P.mean_movement, 'Color', [0.6, 0.2,0]);
sdfill(P.time(1:P.single_trial_length), P.mean_movement, std(P.antennal_movement,[],1), [0.6, 0.2,0], "")

patch(x,y_stim,[0 0 0], 'FaceAlpha' , 0.3, 'EdgeColor', 'none');
% ylabel('Position (deg)', 'Rotation',0);
% yyaxis right; plot(time, P.intendedStimulus(1,:), 'r', 'LineWidth', 0.5);
% ylabel('Generated position stimulus (a.u)');
% title(join([replace([P.date P(1).filename], '_','-')], '   '));
xlim([0 Inf]);
ylim([-1.2 1.2]);
ax1.FontSize = tickLabelSize;
ax1.XLabel.FontSize = labelFontSize;
ax1.YLabel.FontSize = labelFontSize;
ax1.XAxis.Visible = 'off';
ax1.FontName = 'Calibri';
ax1.Box = "off";


ax2 = subplot(2,1,2);
plot(P.time(1:P.single_trial_length), mean(P.gcfr,1), 'Color', [0.4660 0.6740 0.1880]);
sdfill(P.time(1:P.single_trial_length), mean(P.gcfr,1), std(P.gcfr,[],1), [0.4660 0.6740 0.1880], ""); hold on;
% patch(x,y_gcfr,[0 0 0], 'FaceAlpha' , 0.3, 'EdgeColor', 'none');
ax2.FontSize = tickLabelSize;
ax2.FontName = 'Calibri';
ax2.Box = "off";
xlim([0 Inf]);
% ylim([-10 150]);
linkaxes([ax1 ax2], 'x');
savefigures(P(1), "stim_resp", gcf, "png", 'D:\Work\Figures for presentation\stair');


% Split and overlapped stimulus and response

figure;

ax1 = subplot(2,1,1); 
plot(time(1:P.single_trial_length/2), stimulus(1:mid_point), "Color",'k', "LineWidth",1); hold on;
plot(time(1:P.single_trial_length/2), flip(stimulus(mid_point+1:end)), "LineWidth",1, "Color",[0.502,0.502,0.502]);
ylabel('Position (deg)', 'FontSize',labelFontSize, 'Rotation',0);
ylim([-1.2 1.2]);
% title(join([num2str(j) replace([P(1).date P(1).filename], '_','-')], '   '));
ax1.XAxis.Visible = 'off';
ax1.FontSize = tickLabelSize;
ax1.FontName = 'Calibri';
ax1.Box = "off";

ax2 = subplot(2,1,2); plot(time(1:P.single_trial_length/2), P.avg_gcfr(1:mid_point), "Color",'k', "LineWidth",1); hold on;
plot(time(1:P.single_trial_length/2), flip(P.avg_gcfr(mid_point+1:end)), "LineWidth",1, "Color",[0.502,0.502,0.502]);
ylabel('Mean firing rate (Hz)',"FontSize",labelFontSize,"Rotation",0);
xlabel('Time (s)',"FontSize",labelFontSize);
ylim([-10 150]);
ax2.FontSize = tickLabelSize;
ax2.FontName = 'Calibri';
ax2.Box = "off";
xlim([0 Inf]);
linkaxes([ax1 ax2], 'x');


%}

% Hysteresis
S = abs(diff(P.intendedStimulus(1,1:P.single_trial_length)));
x = [];
for i=1:length(S)-1
    if (S(i)>0 && S(i+1)==0)
        x = [x i]; %start point of every static position
    end
end

avg_gcfr = mean(P.gcfr, 1);
x = x+3*P.fs;
x = [3*P.fs x];
positions = P.mean_movement(x);
pos_fr = avg_gcfr(x);

y = [x'-0.5*P.fs x'+0.5*P.fs];

ssFR = nan(P.complete_trials, length(y));
size(ssFR)
for irow = 1:P.complete_trials
    for j = 1:length(y)
        ssFR(irow,j) = mean(P.gcfr(irow, y(j,1):y(j,2)));
    end
end


% Hysteresis plot

err = std(ssFR, [],1);
c = gray(length(positions(4:14)));
pos_fr_onecycle = pos_fr(4:14);
figure;%('WindowState','minimized');

% errorbar(positions, pos_fr, err, 'k-'); hold on;
% scatter(positions, pos_fr, 40, "filled", "CData", c, "MarkerEdgeColor","k");

% errorbar(positions(4:14), pos_fr(4:14), err(4:14), 'k-'); hold on;
% fill(positions(4:14), pos_fr(4:14), [0.4660 0.6740 0.1880], 'FaceAlpha',0.2, 'EdgeColor','none'); hold on;
% scatter(positions(4:14), pos_fr(4:14), 30, "filled", "CData", c, "MarkerEdgeColor","k");

errorbar(positions(4:14), pos_fr_onecycle, err(4:14), 'k-'); hold on;
fill(positions(4:14), pos_fr_onecycle, [0.4660 0.6740 0.1880], 'FaceAlpha',0.2, 'EdgeColor','none'); hold on;
scatter(positions(4:14), pos_fr_onecycle, 30, "filled", "CData", c, "MarkerEdgeColor","k");

% polyin = polyshape(positions(4:14), pos_fr(4:14), 'Simplify',1);
polyin = polyshape(positions(4:14), pos_fr_onecycle, 'Simplify',1);
% polyout = convhull(polyin);
% plot(polyout);
hysteresis_area = area(polyin);

colormap(c)
colorbar(gca, "eastoutside","Box","off", "Ticks", [0 1], "TickLabels", {"Start", "end"}, "Direction","reverse")

xlim([-1.5 1.5]);
ylim([0 Inf]);
yticks(0:25:200);
ylabel('Normalised Mean steady state firing rate (Hz)');
xlabel('Angular position (deg)');
box off;
text(-1,50, string(hysteresis_area));


% savefigures(P(1), "hysteresis_area", gcf, "png", 'D:\Work\Figures for presentation\stair\images');

% Relative position and firing rate
%{
delta_pos = diff(positions);
delta_fr = diff(pos_fr);

figure;
scatter(delta_pos, delta_fr, 'ko', 'filled');
ylabel('Relative change in Firing rate (Hz)');
xlabel('Relative change in position (deg)');
box off;
%}
end