nwb_in = nwbRead('test.nwb');

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


d_rec = designfilt('bandpassiir','FilterOrder',2, ...
'HalfPowerFrequency1',10,'HalfPowerFrequency2',2000, ...
'SampleRate',fs);

d_hes = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',200,'PassbandRipple',0.2, ...
         'SampleRate',10000);

filtered_data_bp = filtfilt(d_rec, rec_data);
hes_data_filtered = filtfilt(d_hes, hes_data);
a= .9258; b=93.15; c=-1.455;
antennal_movement = (b./(hes_data_filtered -a)).^(1/3) + c;
% figure();
% subplot(2,1,1); plot(time, data)
% subplot(2,1,2);



figure();
subplot(3,1,1);plot(time, filtered_data_bp*100, 'k');
A1 = gca;
title('Response to antennal movement')
ylabel('Voltage (mV)')

subplot(3,1,2); hold on; plot(time, antennal_movement, 'Color', [0.6, 0.2,0]);
ylabel('Antennal movement (mm)')
A2 = gca;

subplot(3,1,3); hold on; plot(time, stim_fb,'Color', [0.2,0.3,0.49]);
ylabel('Indenter feedback (V)')
xlabel('time (s)')
A3 = gca;

linkaxes([A1 A2 A3], 'x');

net_movement_antenna = max(antennal_movement) - min(antennal_movement)

%%
rec_protocols = reshape(filtered_data_bp, [no_of_protocols, single_protocol_length]);

rec_chir = filtered_data_bp(1:single_protocol_length);
rec_chirp = reshape(rec_chir, [no_of_trials, single_trial_length]);

rec_amp_sweep = rec_protocols(2,:);
rec_white_noise = rec_protocols(3,:);

% p = [];
% l=[];
% for i = 1:no_of_trials
%     [p(i,:),l(i,:)] = findpeaks(rec_trials(i,:), "MinPeakHeight", 0.02);
% end

%[p,l] = findpeaks(filtered_data_bp, "MinPeakHeight", 0.02);
%l_spikes = zeros(1,length(filtered_data_bp));
%l_spikes(l) = 1;

%l_timestamps = l/fs;

%spiketime_trials = reshape(l_spikes, [3,length(l_spikes)/3]);

figure();
plot(time(1:single_trial_length), 2+3*stim_fb(1:single_trial_length), 'Color', [0.2,0.3,0.49] );
hold on;
% chirp_freq = zeros(1, single_trial_length);
% chirp_freq(3*fs+1:7*fs) = linspace(0,1,4*fs);
% plot(time(1:single_trial_length), chirp_freq, 'Color', 'r');

%rasterplot(l,3,100001,gca , 10000);
raster_data = zeros(no_of_trials, single_trial_length);

k=4;
for i=1:no_of_trials
    [p,l] = findpeaks(rec_chirp(i,:), "MinPeakHeight",0.5*max(rec_chirp(i,:)));
    raster_data(i,l) = 1;
    spike_time = l/fs;
    for j = 1:length(spike_time)
        line([spike_time(j) spike_time(j)], [k k+0.5], 'Color', 'k');
    end
    %plot(spiketime_trials(i,:), '.');
    k = k+1;
end

ylim([0 10]);

sum_of_spikes = sum(raster_data, 1);

L = 500;
alpha = 8;
gauss_win = gausswin(L, alpha);
gcfr = filter(gauss_win, no_of_trials, sum_of_spikes);
plot(time(1:single_trial_length), 9+gcfr/no_of_trials);

figure;
stim_fb_norm = stim_fb/ max(stim_fb);
spike_phase = asin(stim_fb_norm(l));
plot(time(l), rad2deg(spike_phase), '.');
xlabel('time (s)')
ylabel('phase(deg)');

figure;
histogram(rad2deg(spike_phase), 40);
xlabel('Phase(deg)');
ylabel('Frequency');


