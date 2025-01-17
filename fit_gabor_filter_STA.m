%% Fit gabor filters to STA
STA = cell2mat(T_STA.STA);
nrows = height(T_STA);
sigma = nan(nrows,3);
freq = nan(nrows,3);
phase = nan(nrows, 3);
b = nan(nrows,3);
k = nan(nrows,3);
chi2 = nan(nrows,3);
rsquare = nan(nrows, 3);


t = linspace(-100,0,size(STA,2));


gabor_cos = fittype('exp(-0.5 * (x .^ 2 / sigma ^ 2)) .* (b + k *cos(2 * pi * f * x + phi))', 'independent', 'x', 'coefficients', {'sigma', 'f', 'phi', 'b', 'k'});
gabor_sin = fittype('exp(-0.5 * (x .^ 2 / sigma ^ 2)) .* (b + k *sin(2 * pi * f * x + phi))', 'independent', 'x', 'coefficients', {'sigma', 'f', 'phi', 'b', 'k'});
gabor_sin_sin = fittype('(-b*(sigma^-2)*x.*exp(-0.5 * (x .^ 2 / sigma ^ 2))) - (k*(sigma^-2)*x.*exp(-0.5 * (x .^ 2 / sigma ^ 2)).*sin(2*pi*f*x + phi)) + (2*pi*f*k*exp(-0.5 * (x .^ 2 / sigma ^ 2)).* cos(2*pi*f*x + phi))',...
                        'independent', 'x', 'coefficients', {'sigma', 'f', 'phi', 'b', 'k'});
    

gab_cos = @(beta, x) (exp(-0.5 * (x.^2 / beta(1) ^ 2)) .* (beta(4) + beta(5) *cos(2 * pi * beta(2) * x + beta(3))));
gab_sin = @(beta, x) (exp(-0.5 * (x.^2 ./ beta(1) ^ 2)) .* (beta(4) + beta(5) *sin(2 * pi * beta(2) * x + beta(3))));
gab_sin_sin = @(beta, x) ((-beta(4)*(beta(1)^-2)*x.*exp(-0.5 * (x.^2 / beta(1) ^ 2))) - (beta(5)*(beta(1)^-2)*x.*exp(-0.5 * (x .^ 2 / beta(1) ^ 2)).*sin(2*pi*beta(2)*x + beta(3))) + (2*pi*beta(2)*beta(5)*exp(-0.5 * (x .^ 2 / beta(1) ^ 2)).* cos(2*pi*beta(2)*x + beta(3))));





for i = 1:nrows
    
    close all;
    [beta_cos] = improveStartingPts(gab_cos, t', STA(i,:)');
    [beta_sin] = improveStartingPts(gab_sin,  t', STA(i,:)');
    [beta_sin_sin] = improveStartingPts(gab_sin_sin,  t', STA(i,:)');

    fo_cos = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Robust', 'LAR', 'Algorithm', 'Trust-Region', 'MaxIter', 10000, 'StartPoint', beta_cos);
    fo_sin = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Robust', 'LAR', 'Algorithm', 'Trust-Region', 'MaxIter', 10000, 'StartPoint', beta_sin);
    fo_sin_sin = fitoptions('Method', 'NonlinearLeastSquares', ...
    'Robust', 'LAR', 'Algorithm', 'Trust-Region', 'MaxIter', 10000, 'StartPoint', beta_sin_sin);

    [f_cos,gof_cos] = fit(t', STA(i,:)', gabor_cos, fo_cos)
    [f_sin,gof_sin] = fit(t', STA(i,:)', gabor_sin, fo_sin)
    [f_sin_sin,gof_sin_sin] = fit(t', STA(i,:)', gabor_sin_sin, fo_sin_sin)
    chi2_cos = getchi2(STA(i,end-30:end)', f_cos(t(end-30:end)));
    chi2_sin = getchi2(STA(i,end-30:end)', f_sin(t(end-30:end)));
    chi2_sin_sin = getchi2(STA(i,end-30:end)', f_sin_sin(t(end-30:end)));
    sigma(i,:) = [f_cos.sigma f_sin.sigma f_sin_sin.sigma];
    freq(i,:) = [f_cos.f f_sin.f f_sin_sin.f];
    phase(i,:) =[f_cos.phi f_sin.phi f_sin_sin.phi];
    b(i,:) = [f_cos.b f_sin.b f_sin_sin.b];
    k(i,:) = [f_cos.k f_sin.k f_sin_sin.k];
    chi2(i,:) = [chi2_cos chi2_sin chi2_sin_sin];
    rsquare(i,:) = [gof_cos.rsquare gof_sin.rsquare gof_sin_sin.rsquare];
    % 
    % figure; plot(f_cos, t, STA(i,:)); title('Cosine gabor fit'); subtitle(replace(join([T_STA.date(i) T_STA.filename(i)], " "), "_", " "));
    % figure; plot(f_sin, t, STA(i,:)); title('Sine gabor fit'); subtitle(replace(join([T_STA.date(i) T_STA.filename(i)], " "), "_", " "));
    % figure; plot(f_sin_sin, t, STA(i,:)); title('Sine*Sine gabor fit'); subtitle(replace(join([T_STA.date(i) T_STA.filename(i)], " "), "_", " "));
    % 
    % pause;

end

T_STA.sigma = sigma;
T_STA.freq = freq;
T_STA.phase = phase;
T_STA.chi2 = chi2;
T_STA.mean_shift = b;
T_STA.cos_sin_amp = k;
T_STA.rsquare = rsquare;

%% 
% 
T_fit_params = T_STA(:,1:5);
cellID = replace(join([T_STA.date T_STA.filename], " "),"_"," ");
cellID = extractBetween(cellID, "2022-"," T" | " blwgn");
T_fit_params.cellID = cellID;

for i=1:height(T_STA)
    
    [val,idx] = min(T_STA.chi2(i,:));
    T_fit_params.sigma(i) = T_STA.sigma(i,idx);
    T_fit_params.freq(i) = T_STA.freq(i,idx);
    T_fit_params.phase(i) = T_STA.phase(i,idx);
    T_fit_params.mean_shift(i) = T_STA.mean_shift(i,idx);
    T_fit_params.cos_sin_amp(i) = T_STA.cos_sin_amp(i,idx);


end

%% Clustering STAs based on fit parameters
parameters = table2array(T_fit_params(:,7:11), 'AsArray', 'true');
eva = evalclusters(parameters, 'kmeans', 'gap', 'KList',1:6);
ngrp = eva.OptimalK;
[cidx, ctrs] = kmeans(parameters,ngrp,'dist','corr','rep',10, ...
    'MaxIter', 10000, 'disp','final', 'emptyAction', 'drop');
figure;
colororder(parula);
for c = 1:ngrp
    ax(c) = subplot(ngrp,1,c);
    plot(t(end-500:end),STA((cidx == c),end-500:end)', 'LineWidth',1); hold on;
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

%% PCA
[pc, zscores, pcvars] = pca(parameters);

% [coeff,score,latent,tsquared,explained,mu] = pca(parameters);

figure;
stairs(cumsum(pcvars./sum(pcvars)));

figure
pcclusters = clusterdata(zscores(:,1:end),'maxclust',3,'linkage','av');
gscatter(zscores(:,1),zscores(:,2),pcclusters,parula(3));
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');

figure;
scatter3(zscores(:,1), zscores(:,2), zscores(:,3), 10,'filled');
%%

function chi2 = getchi2(O,E)
    chi2  = sum(((O-E).^2)./E);
end

function beta_refined = improveStartingPts(model, xData, yData)
    
    b=1; sigma = 1; f=0.1; phi = 0; k=1;
    options = optimset('MaxFunEvals', 10000);
    beta = [sigma, f, phi, b, k];
    beta_refined = fminsearch(@(beta, x) norm(yData - model(beta, xData)), beta, options);


end
