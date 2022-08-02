function [ev1, ev2] = cov_analysis(raster, antennal_movement, window, fs,start_stim,stop_stim)

    % function [cov_matrix, ev1, ev2] = cov_analysis(raster_data, stimulus, window, fs)
% 
raster_data = raster(:,start_stim:stop_stim);
stimulus = antennal_movement(:,start_stim:stop_stim);
% window = 0.04;


[m,~] = size(raster_data);
STA_freq = [];
all_spike_triggers = [];
for i=1:m
    spike_locs = find(raster_data(i,:)==1);
    if isempty(spike_locs)
        continue;
    end
    for j=1:length(spike_locs)
        if (spike_locs(j)-window*fs)>= 0
            
            spike_triggers = stimulus(1,(spike_locs(j)-window*fs):spike_locs(j));
            all_spike_triggers = [all_spike_triggers; spike_triggers];
        end
    end

end



STA = mean(all_spike_triggers,1);
% t= linspace(-window*1000, 0, length(STA));
% figure;
% plot(t, STA);
% title('Spike Triggered Average');
% xlabel('Time before spike (ms)');
% ylabel('Antennal movement');
% ax = gca;
% ax.Box = "off";


% cov_matrix = cov((all_spike_triggers-STA)');
nspikes = size(all_spike_triggers,1);
cov_matrix = (1/(nspikes -1))*(all_spike_triggers-STA)'*(all_spike_triggers-STA);

% figure;
% 
% h = heatmap((cov_matrix), 'Colormap', parula);
% h.YDisplayData = flipud(h.YDisplayData);
% h.GridVisible = "off";
% Labels = linspace(-window*1e3,0,length(cov_matrix));
% CustomLabels = string(Labels);
% CustomLabels(mod(Labels,10) ~= 0) = " ";
% h.XDisplayLabels = CustomLabels;
% h.YDisplayLabels = flip(CustomLabels);
% title("Spike Triggered Covariance (STC)");

pattern_length = size(all_spike_triggers,2);
m = size(all_spike_triggers,1);
stimulus_prior = zeros([m pattern_length]);
r = randi([1 (length(stimulus)-pattern_length)],1,m);
for j=1:m
    stimulus_prior(j,:) = stimulus(randperm(size(stimulus,1),1),r(j):r(j)+pattern_length-1);
end

avg_stim = mean(stimulus_prior,1);
% stim_prior_cov = cov(stimulus_prior);
stim_prior_cov = (1/(m -1))*(stimulus_prior-avg_stim)'*(stimulus_prior-avg_stim);
% figure;
% h = heatmap((stim_prior_cov), 'Colormap', parula);
% h.YDisplayData = flipud(h.YDisplayData);
% h.GridVisible = "off";
% Labels = linspace(-window*1e3,0,length(cov_matrix));
% CustomLabels = string(Labels);
% CustomLabels(mod(Labels,10) ~= 0) = " ";
% h.XDisplayLabels = CustomLabels;
% h.YDisplayLabels = flip(CustomLabels);
% title("Random stimulus prior");


diff_cov = (cov_matrix - stim_prior_cov) ; 
% figure;
% h = heatmap(diff_cov, 'Colormap', parula);
% h.YDisplayData = flipud(h.YDisplayData);
% h.GridVisible = "off";
% Labels = linspace(-window*1e3,0,length(cov_matrix));
% CustomLabels = string(Labels);
% CustomLabels(mod(Labels,10) ~= 0) = " ";
% h.XDisplayLabels = CustomLabels;
% h.YDisplayLabels = flip(CustomLabels);
% title("STC-C");


[V,D,W] = eig(diff_cov);

[~, ind] = sort(diag(D), 'descend');
Ds = D(ind,ind);
Vs = V(:,ind);
Ws = W(:,ind);

ev1 = Vs(:,1);
ev2 = Vs(:,2);

% figure;
% plot(diag(Ds), '.')
% title('Eigen values');
% ylabel('Eigen values');
% xlabel('Index');
% ax = gca;
% ax.Box = "off";


t= linspace(-window*1000, 0, length(STA));

diag_Ds = diag(Ds);
% plot(t, Vs(:,1).* diag_Ds(1), t, Vs(:,2).*diag_Ds(2));
% ax = gca;
% ax.Box = 'off';
% title('Eigen vectors of STC-C');
% xlabel('Time before spike (ms)');
% ylabel('Antennal movement');




% [U,S,W] = svd(diff_cov, 'econ');
% D_svd = diag(S);
% 
% figure;
% % plot(t, W(:,1:2));
% plot(t, W(:,1).*D_svd, t, W(:,2).*D_svd);
% % title('Eigen vectors of diff cov by SVD');
% ax = gca;
% ax.Box = 'off';
% xlabel('Time before spike (ms)');
% title ('Eigen vectors by SVD');
% ylabel('Antennal movement');

% figure;
% plot(t, U(:,1).*D(1).*W(1,:)', t, U(:,2).*D(2).*W(2,:)');

end




