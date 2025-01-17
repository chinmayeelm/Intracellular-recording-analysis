%% Prediction with NLD
%% Obtain spike triggering ensembles

% linearFilter = sig_evec(:,1);

stim_window = 0.03;
stimulus = P(1).antennal_movement(:,startPt:stopPt);
raster_data = P(1).raster(:,startPt:stopPt);
fs = P(1).fs;
ON_dur = P(1).ON_dur;
linearFilter = STA;%(end-stim_window*fs:end);
% linearFilter = ev(end-0.02*fs:end,2);

numSpikes = length(find(raster_data));
nTrials = P(1).complete_trials;
stimulus_test = stimulus(:,8*fs+1:end);
raster_data_model = raster_data(:,1:8*fs);
raster_data_test = raster_data(:,8*fs+1:end);
gcfr_test = P.gcfr(:,8*fs+1:end);
mean_gcfr = P.avg_gcfr(startPt:stopPt);
mean_gcfr_model = mean_gcfr(1:8*fs);
gcfr = P(1).gcfr(:,startPt:stopPt);
gcfr_test_set = gcfr(1,8*fs+1:end);
mean_FR = mean(mean_gcfr_model); %numSpikes/(nTrials*ON_dur);

spike_window = P.fs*1e-3;

locs_ref = find(raster_data_model(1,:)==1); %spikes of first trial only
ref = locs_ref(1);

stim_pattern = zeros(length(locs_ref),length(linearFilter));
Pc_spike_stim = nan(1,length(locs_ref));
FR_spike_stim = nan(1,length(locs_ref));
similarity_STA = nan(length(locs_ref),1);
% similarity_EVec1 = nan(length(locs_ref),1);
% similarity_EVec2 = nan(length(locs_ref),1);



% Nspikes = length(locs_ref);

% time_train_set = time(1:8*fs);
% figure;
% ax1 = subplot(3,1,1); plot(time_train_set, mean(stimulus_model,1)); hold on;
% ylabel('Avg stimulus (deg)');
% ax2 = subplot(3,1,2); hold on;
% ylabel('Trials #');
% ax3 = subplot(3,1,3); plot(time_train_set, mean_gcfr_model); hold on;
% ylabel('Firing rate (Hz)');
% xlabel('Time (s)');
% linkaxes([ax1,ax2,ax3], 'x');

nTrials = size(raster_data,1);


numSpikes = sum(raster_data_test,'all');
all_spike_triggers = nan([numSpikes stim_window*fs]);
k=1;


for i=1:nTrials
    spike_locs = find(raster_data_test(i,:)==1);
    if isempty(spike_locs)
        continue;
    end
    for j=1:length(spike_locs)
        if (spike_locs(j)-stim_window* fs)> 0
            
            spike_triggers = stimulus_test(i,(spike_locs(j)-stim_window* fs+1):spike_locs(j));

            all_spike_triggers(k,:) = spike_triggers;
            k=k+1;
        end
    end

end
TF=isnan(all_spike_triggers(:,1));
all_spike_triggers(TF,:)=[];


% Prior covariance matrix
pattern_length = size(all_spike_triggers,2);
NstimPriors = size(all_spike_triggers,1)*100;
stimulus_prior = zeros([NstimPriors pattern_length]);
% [rows, cols] = find(raster_data_test==0);
% idx = cols <= pattern_length;
% rows(idx) = [];
% cols(idx) = [];
% 
% r = randi(length(rows),1,NstimPriors);
% r_row = rows(r);
% r_col = cols(r);
% 
% parfor j=1:NstimPriors
%     stimulus_prior(j,:) = stimulus(r_row(j),r_col(j)-pattern_length+1:r_col(j));
% end
%%
rr = randi([1 (length(stimulus_test)-pattern_length)],1,NstimPriors);
parfor j=1:NstimPriors
    stimulus_prior(j,:) = stimulus_test(randperm(size(stimulus_test,1),1),rr(j):rr(j)+pattern_length-1);
end
     
%%
for j = 1:size(all_spike_triggers,1)
            
            similarity_STA(j) = dot(linearFilter,(all_spike_triggers(j,:))')/(norm(all_spike_triggers(j,:))*norm(linearFilter));%/dot(linearFilter,stimulus_model_prior(1,:));
%             similarity_EVec1(j) = dot(eVec1Norm, (stim_pattern(j,:))');
%             similarity_EVec2(j) = dot(eVec2Norm, (stim_pattern(j,:))');
            % similarity_PSE(j) = dot(linearFilter,(stimulus_prior(j,:))')/(norm(stimulus_prior(j,:))*norm(linearFilter));
        
end
similarity_PSE = nan(size(stimulus_prior));
for j = 1:size(stimulus_prior,1)
    similarity_PSE(j) = dot(linearFilter,(stimulus_prior(j,:))')/(norm(stimulus_prior(j,:))*norm(linearFilter));
end

%%


edges = -1:0.005:1;

figure;
histogram(similarity_STA,'BinEdges',edges);
xlabel('Similarity to linear filter');
ylabel('Counts');
title('Normalized STE histogram');

figure;
histogram(similarity_PSE,'BinEdges',edges);
xlabel('Similarity to linear filter');
ylabel('Counts');
title('Normalized PSE histogram');

[N_STE,ed] = histcounts(similarity_STA, edges); %, 'Normalization','pdf'); 
[N_PSE,~] = histcounts(similarity_PSE, edges); %, 'Normalization','pdf');

N_STE_norm = N_STE/sum(N_STE);
N_PSE_norm = N_PSE/sum(N_PSE);

N = N_STE_norm./N_PSE_norm; %/numSpikes
BinCenters = edges(1:end-1)+(diff(ed)/2);
ind = find(N_STE < 10 | N_PSE < 10 | isnan(N));
N(ind) = [];
BinCenters(ind) = [];
FR = N*mean_FR;

figure;
bar(BinCenters,N);
%% 
%% NLD function
figure;
beta0 = [max(FR), 1, 1];
ft = 'a/(1+exp(-b*(x-c)))';
model = @(beta, x) beta(1)/(1+exp(-beta(2)*(x-beta(3)))); 
options = optimset('MaxFunEvals', 10000);
beta0_refined = fminsearch(@(beta) norm(FR - model(beta, BinCenters')), beta0, options);
% beta0_refined = fminsearch(@(beta) norm(FR - fun(beta, BinCenters)), beta0);

[f_pspike, gof] = fit(BinCenters', FR',ft, 'start', beta0_refined) %, 'Exclude', FR<0.001 | isnan(FR));
plot(f_pspike,BinCenters,FR','ko');
xlabel('Projection onto linear filter');
ylabel('Mean firing rate (spikes/s)');
legend('Location', 'best');
title(replace([P(1).date P(1).filename], '_','-'));
hold on;
text(0.2,400,string(["rsquare =" gof.rsquare]));
gof


%% Prediction with projection

predStepWin = 1;%5e-3*fs; %1e-3*fs; % 1 ms

time_test_set = linspace(1,2+stim_window,(2+stim_window)*fs);
stimulus_test_set = [zeros(1,stim_window*fs) stimulus_test(1,:)]; 
figure;
ax1 = subplot(2,1,1); plot(time_test_set, stimulus_test_set); hold on;
ylabel('Avg stimulus (deg)');

pspikePred = [];
timestamps = [];

for iwin = 1:predStepWin:(length(stimulus_test_set)-stim_window*fs)
    test_stim = (stimulus_test_set(1,iwin:iwin+stim_window*fs));
    % subplot(2,1,1); plot(time_test_set((iwin:iwin+stim_window*fs)), test_stim); hold on;
    similarity = dot(test_stim, linearFilter)/(norm(test_stim)*norm(linearFilter));
    timestamps = [timestamps iwin+stim_window*fs];
    pspikePred = [pspikePred f_pspike(similarity)];
end
box off;
set(gca().XAxis, 'Visible', 'off');
% pspikePred = sgolayfilt(pspikePred, 3, (fs/50)+1);
ax2 = subplot(2,1,2);plot(time_test_set(timestamps), pspikePred);%, 'Marker',".","MarkerEdgeColor",'r' );
ylabel('Firing rate');
hold on;
plot(time_test_set(stim_window*fs+1:end),gcfr_test_set, 'Color','#D95319','LineWidth',1);
% ylabel('Actual Firing rate');
title('Predicted firing rate with projection on linear filter')
box off;
legend('Predicted', 'Actual');
legend('Location','best', 'Box','off')
xlabel('Time (s)');

pspikePred_cut = pspikePred(stim_window*fs+1:end);
corrcoef(pspikePred_cut, gcfr_test_set)

%% Actual spike probability for test dataset

% locs_ref = find(raster_data_test(1,:)==1); %spikes of first trial only
% ref = locs_ref(1);
% FR_actual = [];
% 
% for j = 1:length(locs_ref)
%         ref = locs_ref(j);
% 
%         if ref>stim_window*P.fs && (ref+spike_window)<length(stimulus_model)
%             [spike_rows_trials, locs] = find(raster_data_test(:,ref-spike_window : ref+spike_window)==1);
%             FR_actual(j) = mean_gcfr_test_set(j);
% 
%         end      
% end
% % subplot(2,1,2); plot(locs_ref/fs, FR_actual, 'Marker',".","MarkerEdgeColor",'k' );
% ax2 = subplot(2,1,2); plot(time_test_set, gcfr_test, 'k');
% linkaxes([ax1,ax2], 'x');
% legend('Predicted','Observed');