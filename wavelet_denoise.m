% Load the noisy signal
noisy_signal = unfilt;

% Parameters for wavelet denoising
wname = 'db4'; % Wavelet name (choose according to your preference)
level = 2; % Decomposition level
thresholding_method = 'sqtwolog'; % Thresholding method

% Perform wavelet decomposition
[c, l] = wavedec(noisy_signal, level, wname);

% Estimate noise level using the reference baseline noise (assuming it's available)
noise_level = std(baseline_unfilt);

% Apply soft thresholding to coefficients
for i = 1:length(c)
    % Apply thresholding only to detail coefficients (not approximation)
    if isequal(mod(i,level+1),0)
        continue;
    end
    % Soft thresholding
    c(i) = wthresh(c(i), 's', noise_level*sqrt(2*log(length(c))));
end



% Reconstruct the denoised signal
denoised_signal = waverec(c, l, wname);

% Plot the original and denoised signals
figure;
subplot(2,1,1);
plot(noisy_signal);
title('Noisy Signal');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(denoised_signal);
title('Denoised Signal');
xlabel('Time');
ylabel('Amplitude');
