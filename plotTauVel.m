function plotTauVel(T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
uniqueGroups = unique(T.neuronID);
c = parula(length(uniqueGroups));
figure;
hold on;
for i = 1:length(uniqueGroups)
    groupData = T(T.neuronID == uniqueGroups(i), :);
    plot(groupData.velocity, abs(groupData.Tphasic), 'Color', c(i,:), 'Marker','o', 'DisplayName', string(groupData.neuronID(1)));
    
end
hold off;
ylabel('Phasic {\tau} (s)');
xlabel('Velocity ({\circ}/s)')
ax=gca;
ax.YAxis.Scale = 'log';

figure;
hold on;
for i = 1:length(uniqueGroups)
    groupData = T(T.neuronID == uniqueGroups(i), :);
    plot(groupData.velocity, abs(groupData.slopePhasic), 'Color', c(i,:), 'Marker','o', 'DisplayName', string(groupData.neuronID(1)));
    
end
hold off;
ylabel('Rate of change of Firing rate (Hz/s)');
xlabel('Velocity ({\circ}/s)')
ax=gca;
ax.YAxis.Scale = 'linear';
end