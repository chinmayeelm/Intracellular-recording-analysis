function [f_pspike_STA, gof] = constructNLD(STE, PSE, linearFilter, mean_FR)
%CONSTRUCTNLD This function returns the nonlinear decision function
%   Detailed explanation goes here

nFilters = size(linearFilter,1);
similarity_STE = nan(nFilters,size(STE,1));

for i = 1:nFilters
for j = 1:size(STE,1)
            similarity_STE(i,j) = dot(linearFilter(i,:),(STE(j,:))')/(norm(STE(j,:))*norm(linearFilter(i,:)));
end
end

similarity_PSE = nan(nFilters,size(PSE,1));
for i = 1:nFilters
for j = 1:size(PSE,1)
    similarity_PSE(i,j) = dot(linearFilter(i,:),(PSE(j,:))')/(norm(PSE(j,:))*norm(linearFilter(i,:)));
end
end

binWidth = 0.05*std(similarity_PSE(1,:));
edges = -1:binWidth:1;


figure;
for i =1:nFilters
histogram(similarity_STE(i,:),'BinEdges',edges, 'FaceAlpha',0.2, 'Normalization','probability'); hold on;
end
xlabel('Similarity to linear filter');
ylabel('Counts');
title('Normalized STE and PSE histogram');

% figure;
for i=1:nFilters
histogram(similarity_PSE(i,:),'BinEdges',edges, 'FaceAlpha',0.2, 'Normalization','probability'); hold on;
end
xlabel('Similarity to linear filter');
ylabel('Counts');
% title('Normalized PSE histogram');



[N_STE_STA,ed] = histcounts(similarity_STE(1,:), edges); %, 'Normalization','pdf'); 
[N_PSE_STA,~] = histcounts(similarity_PSE(1,:), edges); %, 'Normalization','pdf');

N_STE_norm_STA = N_STE_STA/sum(N_STE_STA);
N_PSE_norm_STA = N_PSE_STA/sum(N_PSE_STA);

N_STA = N_STE_norm_STA./N_PSE_norm_STA; %/numSpikes
BinCenters_STA = edges(1:end-1)+(diff(ed)/2);
ind = find(N_STE_STA < 20 | N_PSE_STA < 20 | isnan(N_STA));
N_STA(ind) = [];
BinCenters_STA(ind) = [];
FR = N_STA.*mean_FR;

figure;
bar(BinCenters_STA,N_STA);

% NLD function
%{
beta0 = [max(FR), 1, 1];
ft = 'a/(1+exp(-b*(x-c)))';
model = @(beta, x) beta(1)/(1+exp(-beta(2)*(x-beta(3)))); 
options = optimset('MaxFunEvals', 10000);
beta0_refined = fminsearch(@(beta) norm(FR - model(beta, BinCenters')), beta0, options);
[f_pspike, gof] = fit(BinCenters', FR',ft, 'start', beta0_refined) %, 'Exclude', FR<0.001 | isnan(FR));
%}

[f_pspike_STA, gof] = fit(BinCenters_STA', FR', 'poly3');



figure;
plot(f_pspike_STA,BinCenters_STA,FR','ko');
xlabel('Projection onto linear filter');
ylabel('Mean firing rate (spikes/s)');
legend('Location', 'best');
hold on;
text(0.2,400,string(["rsquare =" gof.rsquare]));


if nFilters > 1
% figure;
% histogram2(similarity_STE(end-1,:), similarity_STE(end,:), 'XBinEdges',edges,'YBinEdges',edges);
% xlabel('Linear filter 1')
% ylabel('Linear filter 2');

%{
figure;
histogram2(similarity_PSE(end-1,:), similarity_PSE(end,:), 'XBinEdges',edges,'YBinEdges',edges);
xlabel('Linear filter 1');
ylabel('Linear filter 2');
title('PSE');
%}

[N_STE_EV, Xedges, Yedges] = histcounts2(similarity_STE(end-1,:), similarity_STE(end,:), 'XBinEdges',edges,'YBinEdges',edges);
[N_PSE_EV, ~, ~] = histcounts2(similarity_PSE(end-1,:), similarity_PSE(end,:), 'XBinEdges',edges,'YBinEdges',edges);

N_STE_EV_norm = N_STE_EV/sum(N_STE_EV,'all');
N_PSE_EV_norm = N_PSE_EV/sum(N_PSE_EV,'all');

N_EV = N_STE_EV_norm./N_PSE_EV_norm;

FR_ev = N_EV.*mean_FR;

% figure;
% surf(FR_ev);



end

end