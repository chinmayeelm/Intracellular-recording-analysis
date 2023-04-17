function sdfill(time,meanMovement, sd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

curve1 = meanMovement + sd;
curve2 = meanMovement - sd;
x2 = [time fliplr(time)];
inBetween = [curve1 fliplr(curve2)];
fill(x2, inBetween, 'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
hold on;
plot(time, meanMovement, 'k', 'LineWidth', 0.5);
end

