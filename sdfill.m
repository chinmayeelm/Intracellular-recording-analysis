function sdfill(time,meanMovement, sd, c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

curve1 = meanMovement + sd;
curve2 = meanMovement - sd;
time = reshape(time, 1,[]);
x2 = [time fliplr(time)];
inBetween = [curve1 fliplr(curve2)];
hold on;
fill(x2, inBetween, c, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

plot(time, meanMovement, 'Color', c, 'LineWidth', 1); 
% plot(time, curve1, 'Color', c, 'LineWidth', 1);
% plot(time, curve2, 'Color', c, 'LineWidth', 1);
end

