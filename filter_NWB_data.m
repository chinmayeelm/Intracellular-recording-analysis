nwb_in = nwbRead('M1_N1_T2.nwb');

disp(nwb_in)

rec = nwb_in.acquisition.get('response_to_JO_stimulation');
data = rec.data.load;
rec_data = data(:,1);
stim_fb = data(:,2);
hes_data = data(:,3);
t = rec.timestamps.load;

no_of_protocols = 3;

no_of_trials = 5;
single_protocol_length = length(rec_data)/no_of_protocols;

single_trial_length = single_protocol_length/no_of_trials;

fs = 10000; %sampling freq
time = t;
amplification_fator = 10;
milivolts_conv_factor = 1000/amplification_fator;


d = designfilt('bandpassiir','FilterOrder',2, ...
'HalfPowerFrequency1',10,'HalfPowerFrequency2',2000, ...
'SampleRate',fs);


filtered_data_bp = filtfilt(d, rec_data)*milivolts_conv_factor;

a= .9258; b=93.15; c=-1.455;
antennal_movement = (b./(hes_data -a)).^(1/3) + c;
% figure();
% subplot(2,1,1); plot(time, data)
% subplot(2,1,2);



net_movement_antenna = max(antennal_movement) - min(antennal_movement)
rec_protocols = reshape(filtered_data_bp,[single_protocol_length, no_of_protocols]);
stim_protocols = reshape(antennal_movement, [single_protocol_length, no_of_protocols]);
%%
%rec_freq_chirp = filtered_data_bp(1:single_protocol_length);
rec_freq_chirp = reshape(rec_protocols(:,1), [single_trial_length, no_of_trials]);
chirp_stim = stim_protocols(1:single_trial_length,1);
%rec_chirp = rec_chirp_';

rec_amp_sweep = reshape(rec_protocols(:,2), [single_trial_length, no_of_trials]);
amp_stim = stim_protocols(1:single_trial_length,2);

rec_white_noise = reshape(rec_protocols(:,3), [single_trial_length, no_of_trials]);
gwn_stim = stim_protocols(1:single_trial_length,3);


plot_figures(fs, time, single_trial_length, chirp_stim, rec_freq_chirp,no_of_trials); 
plot_figures(fs, time, single_trial_length, amp_stim, rec_amp_sweep,no_of_trials); 
plot_figures(fs, time, single_trial_length, gwn_stim, rec_white_noise,no_of_trials); 

function f = plot_figures(fs, time, single_trial_length, antennal_movement, rec_trace, no_of_trials)

    figure();
    subplot(4,1,1); hold on; plot(time(1:single_trial_length), antennal_movement(1:single_trial_length), 'Color', '#D95319' );
    ylabel('Antennal movement (mm)');
    %ylim([0.9*min(antennal_movement(1:single_trial_length)) 1.1*max(antennal_movement(1:single_trial_length))]);
    set(gca, 'xtick', []);
    A1 = gca;


    [pks,locs] = findpeaks(rec_trace(:,1), "MinPeakHeight",0.5*max(rec_trace(:,1)));
    spike_locs = locs/fs;
    subplot(4,1,2);plot(time(1:single_trial_length), rec_trace(:,1), spike_locs, pks, '.');
    ylim([1.25*min(rec_trace(:,1)) 1.25*max(pks)]);
    set(gca, 'xtick', []);
    A2 = gca;
    % title('Response to antennal movement')
    ylabel('Voltage (mV)')



    %%
    raster_data = zeros(no_of_trials, single_trial_length);

    k=0.3;
    for i=1:no_of_trials
        [p,l] = findpeaks(rec_trace(:,i), "MinPeakHeight",0.5*max(rec_trace(:,i)));
        raster_data(i,l) = 1;
        spike_time = l/fs;
        for j = 1:length(spike_time)
            subplot(4,1,3); line([spike_time(j) spike_time(j)], [k k+0.5],'Color', 'k');
            ylim([0 5]);
        end
        %plot(spiketime_trials(i,:), '.');
        k = k+1;
    end

    ylabel('Trials');
    set(gca, 'xtick', []);
    A3 = gca;

    sum_of_spikes = sum(raster_data, 1);

    L = 500;
    alpha = 8;
    gauss_win = gausswin(L, alpha);
    gcfr = filter(gauss_win, no_of_trials, sum_of_spikes);
    subplot(4,1,4); plot(time(1:single_trial_length), gcfr/no_of_trials, 'Color', '#EDB120' );
    ylim([-0.1 1.1*max(gcfr/no_of_trials)]);
    A4 = gca;

    ylabel('Avg. GCFR');
    xlabel('time (s)');

    linkaxes([A1 A2 A3 A4], 'x');

end
% figure;
% stim_fb_norm = stim_fb/ max(stim_fb);
% spike_phase = asin(stim_fb_norm(l));
% plot(time(l), rad2deg(spike_phase), '.');
% xlabel('time (s)')
% ylabel('phase(deg)');
% 
% figure;
% histogram(rad2deg(spike_phase), 40);
% xlabel('Phase(deg)');
% ylabel('Frequency');

 
