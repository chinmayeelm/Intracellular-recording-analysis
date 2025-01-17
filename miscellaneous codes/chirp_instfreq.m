% Define parameters
fs = P(end).fs;          % Sampling frequency (Hz)
t = linspace(0,15,length(stim));       % Time vector (5 seconds)
% f0 = 0;            % Starting frequency of chirp (Hz)
% f1 = 150;           % Ending frequency of chirp (Hz)
% A = 1;              % Amplitude
% SNR = 10;           % Signal-to-Noise Ratio (dB)

% Generate chirp signal
% x_clean = chirp(t, f0, t(end), f1);
% x_noisy = awgn(x_clean, SNR, 'measured');
x_clean = stim_ifb;
x_noisy = stim;

% Compute Hilbert transform to get instantaneous frequency
hilbert_transform = hilbert(x_noisy);
inst_phase = unwrap(angle(hilbert_transform)); % Unwrap phase
inst_freq = diff(inst_phase) * fs / (2 * pi); % Instantaneous frequency

% Plot results
figure;
subplot(3,1,1);
plot(t, x_clean, 'b', t, x_noisy, 'r');
title('Noisy Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Clean Signal', 'Noisy Signal');

subplot(3,1,2);
plot(t, inst_phase);
title('Instantaneous Phase');
xlabel('Time (s)');
ylabel('Phase (rad)');

subplot(3,1,3);
plot(t(1:end-1), inst_freq);
title('Instantaneous Frequency');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 200]);

% Peak detection for better visualization
[~,locs] = findpeaks(inst_freq);
hold on;
plot(t(locs), inst_freq(locs), 'ro');
hold off;


%% Using CWT

% Define parameters
fs = 1000;          % Sampling frequency (Hz)
t = 0:1/fs:5;       % Time vector (5 seconds)
f0 = 10;            % Starting frequency of chirp (Hz)
f1 = 190;           % Ending frequency of chirp (Hz)
A = 1;              % Amplitude
SNR = 10;           % Signal-to-Noise Ratio (dB)

% Generate chirp signal
x_clean = chirp(t, f0, t(end), f1);
x_noisy = awgn(x_clean, SNR, 'measured');

% Perform Continuous Wavelet Transform (CWT)
scales = 1:128;
cwt_signal = cwt(x_noisy, scales, 'cmor1-1');

% Compute instantaneous frequency from CWT
delta_t = 1/fs;
time = (0:length(x_noisy)-1) * delta_t;
inst_freq = scal2frq(cwt_signal, 'cmor1-1', delta_t);

% Plot results
figure;
subplot(3,1,1);
plot(t, x_clean, 'b', t, x_noisy, 'r');
title('Noisy Chirp Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Clean Signal', 'Noisy Signal');

subplot(3,1,2);
imagesc(time, scales, abs(cwt_signal));
title('Continuous Wavelet Transform (CWT)');
xlabel('Time (s)');
ylabel('Scale');
colorbar;

subplot(3,1,3);
imagesc(time, scales, inst_freq);
title('Instantaneous Frequency');
xlabel('Time (s)');
ylabel('Scale');
colorbar;


%% Based on midcross

stim_smooth = sgolayfilt(stim,2,101);
[c, midlev] = midcross(stim_smooth, fs, 'Tolerance',5);
t = linspace(0,15, length(stim));
figure;
plot(t, stim_smooth); hold on;
plot(c',midlev,'rx')

c = [0 c];
T_half = diff(c);
freq = 1./(2*T_half);
t_f = T_half./2;
t_f = c(1:end-1)+T_half./2;

% figure; plot(t,stim_smooth); xline(t_f, 'k--');

f = fit(t_f',freq','poly1');
figure; plot(f,t_f,freq);