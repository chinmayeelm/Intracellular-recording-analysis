function figHandle = stairGCFRplots(P)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
total_trial_dur = P.single_trial_length/P.fs;
time = linspace(0,total_trial_dur,P(1).single_trial_length);

stimulus = -P.mean_movement;
mid_point = length(stimulus)/2;
labelFontSize = 14;
tickLabelSize = 12;

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

y_stim = repmat([-1.2 -1.2 1.2 1.2]', 1,length(x));
y_gcfr = repmat([0 0 150 150]', 1,length(x));

figHandle = figure('Color', 'w', 'Visible','off');

% ax1 = subplot(2,1,1); plot(time, stimulus, 'k'); hold on;
% patch(x,y_stim,[0.8500 0.3250 0.0980], 'FaceAlpha' , 0.3, 'EdgeColor', 'none');
% % ylabel('Position (deg)', 'Rotation',0);
% % yyaxis right; plot(time, P.intendedStimulus(1,:), 'r', 'LineWidth', 0.5);
% % ylabel('Generated position stimulus (a.u)');
% % title(join([replace([P.date P(1).filename], '_','-')], '   '));
% xlim([0 Inf]);
% ax1.FontSize = tickLabelSize;
% ax1.XLabel.FontSize = labelFontSize;
% ax1.YLabel.FontSize = labelFontSize;
% ax1.XAxis.Visible = 'off';
% ax1.FontName = 'Calibri';
% ax1.Box = "off";


ax2 = subplot(2,1,2);
[lineOut, ~] = stdshade(P.gcfr,0.2,'k',time); %10 = (fs/L)*gcfr Hz
lineOut.LineWidth  = 0.05;
% patch(x,y_gcfr, [0.8500 0.3250 0.0980], 'FaceAlpha' , 0.3, 'EdgeColor', 'none');
% ylabel('Mean firing rate (Hz)','Rotation',0, 'FontSize',labelFontSize);
% xlabel('Time (s)','FontSize',labelFontSize);
ax2.FontSize = tickLabelSize;
ax2.FontName = 'Calibri';
ax2.Box = "off";
xlim([0 Inf]);
ylim([0 150]);
% linkaxes([ax1 ax2], 'x');

%{
 figure;
    
 ax1 = subplot(2,1,1); plot(time(1:P.single_trial_length/2), stimulus(1:mid_point), "Color",'k', "LineWidth",1); hold on;
 plot(time(1:P.single_trial_length/2), fliplr(stimulus(mid_point+1:end)), "LineWidth",1, "Color",[0.502,0.502,0.502]);
 ylabel('Position (deg)', 'FontSize',labelFontSize, 'Rotation',0);
 % title(join([num2str(j) replace([P(1).date P(1).filename], '_','-')], '   '));
ax1.XAxis.Visible = 'off';
ax1.FontSize = tickLabelSize;
ax1.FontName = 'Calibri';

 ax2 = subplot(2,1,2); plot(time(1:P.single_trial_length/2), P.avg_gcfr(1:mid_point), "Color",'k', "LineWidth",1); hold on;
 plot(time(1:P.single_trial_length/2), fliplr(P.avg_gcfr(mid_point+1:end)), "LineWidth",1, "Color",[0.502,0.502,0.502]);
 ylabel('Mean firing rate (Hz)',"FontSize",labelFontSize,"Rotation",0);
 xlabel('Time (s)',"FontSize",labelFontSize);
ax2.FontSize = tickLabelSize;
ax2.FontName = 'Calibri';
ax2.Box = "off";
xlim([0 Inf]);
linkaxes([ax1 ax2], 'x');
%}
end

