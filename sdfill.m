function sdfill(time,meanMovement, sd, c, varargin)
% sdfill(time,meanMovement, sd, c)
%   Detailed explanation goes here
if nargin > 4
    fill_flag = varargin{end}
else
    fill_flag = "";
end

curve1 = meanMovement + sd;
curve2 = meanMovement - sd;
time = reshape(time, 1,[]);
x2 = [time fliplr(time)];
inBetween = [curve1 fliplr(curve2)];

plot(time, meanMovement, 'Color', c, 'LineWidth', 1); hold on;

if fill_flag == "fill"
    fill(x2, inBetween, c, 'FaceAlpha', 0.4, 'EdgeColor','none');
elseif fill_flag == "none"
    
else
    plot(time, curve1, 'Color', c, 'LineWidth', 1);
    plot(time, curve2, 'Color', c, 'LineWidth', 1);
end
end

