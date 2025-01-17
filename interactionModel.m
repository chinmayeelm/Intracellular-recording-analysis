%% For velocity and tau
mdl = stepwiselm(T_phasic,'interactions')
plotSlice(mdl)
figure; plotEffects(mdl)
figure; plotInteraction(mdl,'neuronID','velocity')
figure; plotInteraction(mdl,'velocity','peakFR', 'predictions')
figure; plotInteraction(mdl,'velocity','neuronID','predictions')
figure; plotInteraction(mdl,'neuronID','peakFR','predictions')

%% For jitter
mdl = stepwiselm(T_jitter_interaction,'interactions', 'Verbose',2)
figure; plotSlice(mdl)
figure; plotEffects(mdl)
figure; plotInteraction(mdl,'slope','amplitude')

figure; plotInteraction(mdl,'slope_sd','amplitude_sd')

figure; plotInteraction(mdl,'neuronID','time_to_prev_spike')

figure; plotInteraction(mdl,'amplitude_sd','slope_sd', 'predictions')
figure; plotInteraction(mdl,'amplitude','slope','predictions')
figure; plotInteraction(mdl,'neuronID','time_to_prev_spike','predictions')