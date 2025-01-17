function ax = colorLineByMagnitude(x,y, cmapName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Determine colormap
cmap = colormap(cmapName);
numColors = size(cmap, 1);

% Plot the time series with color representing magnitude
for i = 1:numel(x)-1
    % Calculate color index based on magnitude
    colorIndex = round(interp1(linspace(min(y), max(y), numColors), 1:numColors, y(i)));
    
    % Plot line segment with appropriate color
    % ax = line([x(i), x(i+1)], [y(i), y(i+1)], 'Color', cmap(colorIndex, :), 'LineWidth', 2);
    ax = line([x(i), x(i+1)], [y(i), y(i+1)], 'Color', cmap(colorIndex, :), 'Marker', '.', 'MarkerSize', 10);
    
end

% Add colorbar
c = colorbar;
c.Label.String = 'Magnitude';


end