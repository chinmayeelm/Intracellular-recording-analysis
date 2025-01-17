function [sig_evec_orth, STA] = cov_analysis_old(raster, antennal_movement, stim_window, fs,start_stim,stop_stim)


raster_data = raster(:,start_stim:stop_stim);
stimulus = antennal_movement(:,start_stim:stop_stim);
% window = 0.04;

nTrials = size(raster_data,1);


[m,~] = size(raster_data);
numSpikes = sum(raster_data,'all');
all_spike_triggers = nan([numSpikes stim_window* fs+1]);
k=1;


for i=1:m
    spike_locs = find(raster_data(i,:)==1);
    if isempty(spike_locs)
        continue;
    end
    for j=1:length(spike_locs)
        if (spike_locs(j)-stim_window* fs)>= 0
            
            spike_triggers = stimulus(i,(spike_locs(j)-stim_window* fs):spike_locs(j));

            all_spike_triggers(k,:) = spike_triggers;
            k=k+1;
        end
    end

end
TF=isnan(all_spike_triggers(:,1));
all_spike_triggers(TF,:)=[];


% Prior covariance matrix
pattern_length = size(all_spike_triggers,2);
NstimPriors = size(all_spike_triggers,1);
stimulus_prior = zeros([NstimPriors pattern_length]);
[rows, cols] = find(raster_data==0);
idx = cols <= pattern_length;
rows(idx) = [];
cols(idx) = [];

r = randi(length(rows),1,NstimPriors);
r_row = rows(r);
r_col = cols(r);

parfor j=1:NstimPriors
    stimulus_prior(j,:) = stimulus(r_row(j),r_col(j)-pattern_length+1:r_col(j));
end

avg_stim = mean(stimulus_prior,1);
stim_prior_cov = cov(stimulus_prior-avg_stim);
% stim_prior_cov = (1/(NstimPriors -1))*(stimulus_prior-avg_stim)'*(stimulus_prior-avg_stim);
% figure;
% h = heatmap((stim_prior_cov), 'Colormap', parula);
% h.YDisplayData = flipud(h.YDisplayData);
% h.GridVisible = "off";
% Labels = linspace(-stim_window*1e3,0,length(stim_prior_cov));
% CustomLabels = string(Labels);
% CustomLabels(mod(Labels,10) ~= 0) = " ";
% h.XDisplayLabels = CustomLabels;
% h.YDisplayLabels = flip(CustomLabels);
% title("Random stimulus prior");

% Calculating STA
% STA_ = mean(all_spike_triggers,1);
STA = mean(all_spike_triggers-avg_stim,1); %% Subtracting mean of all stimulus patterns
t= linspace(-stim_window*1000, 0, length(STA));

    


nspikes = size(all_spike_triggers,1);
cov_matrix = cov(all_spike_triggers-STA);
% cov_matrix = (1/(nspikes -1))*(all_spike_triggers-STA)'*(all_spike_triggers-STA);


% figure;
% 
% h = heatmap((cov_matrix), 'Colormap', parula);
% h.YDisplayData = flipud(h.YDisplayData);
% h.GridVisible = "off";
% Labels = linspace(-stim_window*1e3,0,length(cov_matrix));
% CustomLabels = string(Labels);
% CustomLabels(mod(Labels,10) ~= 0) = " ";
% h.XDisplayLabels = CustomLabels;
% h.YDisplayLabels = flip(CustomLabels);
% title("Spike Triggered Covariance (STC)");
% Null distribution of eigen values
nPatterns = nspikes;
rand_stim = zeros([nPatterns pattern_length]);
Ds_null = [];

for repeats = 1:500
    rr = randi([1 (length(stimulus)-pattern_length)],1,nPatterns);
    parfor j=1:nPatterns
        rand_stim(j,:) = stimulus(randperm(size(stimulus,1),1),rr(j):rr(j)+pattern_length-1);
    end
    
    random_cov = cov(rand_stim-STA);
    
    [~,D] = eig(random_cov-stim_prior_cov);
    % [V,D] = eig(cov_matrix);
    Ds_null= [Ds_null diag(D)];
end

% [~, ind] = sort(real(Ds_null), 'descend');
% [~, ind] = sort(real(diag(D)), 'descend');


% figure;
% plot(sort(real(Ds_null),'ascend'), '.',"MarkerSize",20);
% % plot(sort(real(eVal),'descend'))
% title('Eigen values of null distribution');
% ylabel('Eigen values');
% xlabel('Index');
% ax = gca;
% ax.Box = "off";

min_null = min(Ds_null, [], 'all');
max_null = max(Ds_null, [], 'all');

diff_cov = (cov_matrix - stim_prior_cov) ; 
% figure;
% h = heatmap(diff_cov, 'Colormap', parula);
% h.YDisplayData = flipud(h.YDisplayData);
% h.GridVisible = "off";
% Labels = linspace(-stim_window*1e3,0,length(cov_matrix));
% CustomLabels = string(Labels);
% CustomLabels(mod(Labels,10) ~= 0) = " ";
% h.XDisplayLabels = CustomLabels;
% h.YDisplayLabels = flip(CustomLabels);
% title("STC-C");


[V,D] = eig(diff_cov);
% [V,D] = eig(cov_matrix);

[~, ind] = sort(abs(real(diag(D))), 'descend');
% [~, ind] = sort(real(diag(D)), 'descend');
Ds = D(ind,ind);
Vs = V(:,ind);

eVal = diag(Ds);

% Significant Eigen values
eVal_sig_ind = find(eVal < min_null | eVal > max_null);
eVal_sig = eVal(eVal_sig_ind);

sig_evec = eVal_sig' .* Vs(:,eVal_sig_ind);


% figure;
% plot(sort(real(diag(D)),'ascend'), '.',"MarkerSize",20); hold on;
% yline(min_null, '--k');
% yline(max_null, '--k');
% title('Eigen values of Cs-Cp');
% ylabel('Eigen values');
% xlabel('Index');
% ax = gca;
% ax.Box = "off";


t= linspace(-stim_window*1000, 0, length(STA));

% figure;
% 
% hold on;
% for val_ind=1:length(eVal_sig)
%     plot(t,sig_evec(:,val_ind)*eVal_sig(val_ind), 'LineWidth',1);
% end
% ylabel('Antennal movement (mm)');
% % plot(t,sum(Vs(:,1:2),2));
% % yyaxis right; plot(t,STA); ylabel("STA");
% % legend("EV1","EV2","Sum of EVs", "STA", 'Location',"best");
% % plot(t,Vs(:,1:2));
% ax = gca;
% % ax.YAxis.Visible = 'off'
% ax.Box = 'off';

% title('Significant eigen vectors');
% xlabel('Time before spike (ms)');
% ylabel('Antennal movement (mm)');


%Orthogonalize eigen vectors
sig_evec_orth = [];
for i = 1:length(eVal_sig)
    sig_evec_orth(:,i) = sig_evec(:,i) -  ((dot(sig_evec(:,i), STA)/(norm(STA))^2)*STA');    
end


end



