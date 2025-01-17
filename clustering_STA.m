%% Clustering STAs

idx = T_STA_EV.fs == 2e4;
sta = cell2mat(T_STA_EV.STA(idx));
% sta_norm = miscFuncs.minmaxNorm(sta);
time = cell2mat(T_STA_EV.time(idx));
ngrp = 3;
[sta_idx, C] = kmeans(sta, ngrp, 'Replicates',10, 'MaxIter', 10000); 

figure; hold on;
c = lines(4);

for i = 1:size(sta,1)
    ax(sta_idx(i)) = subplot(ngrp,1,sta_idx(i));
    plot(time(i,:), sta(i,:), 'LineWidth',1, 'Color', c(sta_idx(i),:), ...
        'DisplayName',replace(join([T_STA_EV.date(i) T_STA_EV.filename(i)]," "), "_", "")); 
    hold on;
    box off;
    grid on;
    ax(sta_idx(i)).XAxis.Visible = "off";
    ylim([-0.4 0.4]);
    
end

ax(ngrp).XAxis.Visible = "on";
linkaxes(ax, 'x');

figure;
plot(sta(sta_idx==1,1),sta(sta_idx==1,2),'r.','MarkerSize',12)
hold on
plot(sta(sta_idx==2,1),sta(sta_idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

%% K-means with optimal cluster evaluation
ngrp = 3;
% patterns = isi_all';% miscFuncs.minmaxNorm(sta(:,501:end));
% t = time(end, 501:end);
patterns = miscFuncs.minmaxNorm_Minus1ToPlus1(sta);
colors = parula(size(patterns,1));

eva = evalclusters(patterns, 'kmeans', 'silhouette', 'KList',1:6);
ngrp = eva.OptimalK;
[cidx, ctrs] = kmeans(patterns,ngrp,'dist','corr','rep',10, ...
    'MaxIter', 10000, 'disp','final', 'emptyAction', 'drop');
figure;
colororder(colors);
for c = 1:ngrp
    ax(c) = subplot(ngrp,1,c);
    plot(patterns((cidx == c),:)'); hold on;
    % plot(ctrs(c,:)', 'k', 'LineWidth',1);
    axis padded
    box off;
    % ylim([-0.4 0.4]);
    % yticks(-0.4:0.2:0.4);
end

linkaxes(ax, 'x');
ax(1).XAxis.Visible = "off";
sgtitle('K-Means Clustering of STA');

%% Hierarchical clustering
eva1 = evalclusters(patterns, 'linkage', 'silhouette', 'KList',1:6);
ngrp1 = eva1.OptimalK;
corrDist = pdist(patterns,'corr');
clusterTree = linkage(corrDist,'average');
clusters = cluster(clusterTree,'maxclust',ngrp1);

figure
for c = 1:ngrp1
    subplot(ngrp1,1,c);
    plot(t,patterns((clusters == c),:)');
    axis tight
    box off;

end
sgtitle('Hierarchical Clustering of Profiles');

%% DTW clustering
patterns = miscFuncs.minmaxNorm_Minus1ToPlus1(sta);
maxClust = 3;
lst2clu = {'s','ca1','ca3','ca6'};
S = mdwtcluster(patterns,'maxclust',maxClust,'lst2clu',lst2clu);
IdxCLU = S.IdxCLU;
figure; hold on;
subplot(maxClust,1,1); plot(patterns(IdxCLU(:,1)==1,:)','r')
subplot(maxClust,1,2);plot(patterns(IdxCLU(:,1)==2,:)','g')
subplot(maxClust,1,3);plot(patterns(IdxCLU(:,1)==3,:)','b')
% subplot(maxClust,1,4);plot(patterns(IdxCLU(:,1)==4,:)','m')
hold off
% title('Cluster 1 (Signal) and Cluster 3 (Coefficients)')
equalPART = isequal(IdxCLU(:,1),IdxCLU(:,3))
%% PCA
[pc, zscores, pcvars] = pca(distances);

% pcvars./sum(pcvars) * 100;
% cumsum(pcvars./sum(pcvars) * 100);    

figure;
stairs(cumsum(pcvars./sum(pcvars)));

figure
pcclusters = clusterdata(zscores(:,1:end),'maxclust',6,'linkage','av');
gscatter(zscores(:,1),zscores(:,2),pcclusters,hsv(8));
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');


%% Clustering with Gabor filters

distances = nan([size(STA,1),size(filters,2)]);

for ista = 1:size(STA,1)
    for jfilt = 1:size(filters, 2)
        distances(ista,jfilt) = dtw(STA(ista,:), filters(:,jfilt));
    end
end

% colors = parula(size(STA,1));


eva = evalclusters(distances(:,[1 3 5]), 'kmeans', 'gap', 'KList',1:10);
ngrp = eva.OptimalK;
[cidx, ctrs] = kmeans(distances(:,[1 3 5]),ngrp,'dist','corr','rep',10, ...
    'MaxIter', 10000, 'disp','final', 'emptyAction', 'drop');
figure;
colororder(parula);
for c = 1:ngrp
    ax(c) = subplot(ngrp,1,c);
    plot(STA((cidx == c),:)', 'LineWidth',1); hold on;
    % plot(ctrs(c,:)', 'k', 'LineWidth',1);
    axis padded
    box off;
    ax(c).XAxis.Visible = "off";
    % ylim([-0.4 0.4]);
    % yticks(-0.4:0.2:0.4);
end

linkaxes(ax, 'x');
ax(ngrp).XAxis.Visible = "on";
sgtitle('K-Means Clustering of STA');
