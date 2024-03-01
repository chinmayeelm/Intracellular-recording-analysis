function plotTauPos(T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
uniqueGroups = unique(T.neuronID);
c = parula(length(uniqueGroups));
figure;
hold on;
for i = 1:length(uniqueGroups)
    groupData = T(T.neuronID == uniqueGroups(i), :);
    subgroupData_low = groupData(groupData.position <=0, :);
    subgroupData_high = groupData(groupData.position > 0, :);
    for j=1:height(subgroupData_low)
        plot(subgroupData_low.position, abs(subgroupData_low.Ttonic), 'Color', c(i,:), 'Marker','o', 'DisplayName', string(subgroupData_low.neuronID(1)));
    end
    for j=1:height(subgroupData_high)
        plot(subgroupData_high.position, abs(subgroupData_high.Ttonic), 'Color', c(i,:), 'Marker','o', 'DisplayName', string(subgroupData_high.neuronID(1)));
    end
end
hold off;
ylabel('Tonic {\tau} (s)');
xlabel('Position ({\circ})')
% ax=gca;
% ax.YAxis.Scale = 'log';

% figure;
% hold on;
% for i = 1:length(uniqueGroups)
%     groupData = T(T.neuronID == uniqueGroups(i), :);
% 
%     subgroupData_low = groupData(groupData.position <=0, :);
%     subgroupData_high = groupData(groupData.position > 0, :);
%     for j=1:height(subgroupData_low)
%         plot(subgroupData_low.position, abs(subgroupData_low.slopeTonic), 'Color', c(i,:), 'Marker','o', 'DisplayName', string(subgroupData_low.neuronID(1)));
%     end
%     for j=1:height(subgroupData_high)
%         plot(subgroupData_high.position, abs(subgroupData_high.slopeTonic), 'Color', c(i,:), 'Marker','o', 'DisplayName', string(subgroupData_high.neuronID(1)));
%     end
% end
% hold off;
% ylabel('Rate of change of Firing rate (Hz/s)');
% xlabel('Position ({\circ})')
% ax=gca;
% ax.YAxis.Scale = 'linear';
end