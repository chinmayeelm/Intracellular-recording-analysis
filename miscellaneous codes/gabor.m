% Define parameters
frequency = 0.1;  % Frequency of the Gabor filter
% theta = pi/4;    % Orientation of the Gabor filter (in radians)
sigma = 2;%5;       % Standard deviation of the Gaussian envelope
L = 400;      % Size of the filter

% Create a grid of points
x = linspace(-40, 0, L);
% x = linspace(0, 40, L);

% Compute the cosine and sine Gabor filters

% cosine_gabor = (gabor_filter(x, frequency, sigma,0));
% sine_gabor_shifted = (gabor_filter(x, frequency, sigma, pi/2));  % Phase shifted by pi/2 for sine Gabor
% sine_sine_gabor = (conv(sine_gabor, sine_gabor, 'same'));

cosine_gabor = circshift(gabor_filter(x, frequency, sigma,0), 350,2);
sine_gabor = gabor_filter(x, frequency, sigma, pi/2);
sine_gabor_shifted = circshift(gabor_filter(x, frequency, sigma, pi/2), 350,2);  % Phase shifted by pi/2 for sine Gabor
sine_sine_gabor = circshift(conv(sine_gabor, sine_gabor, 'same'),300, 2);

% Plot the results
figure;
subplot(1, 3, 1);
plot(x, cosine_gabor, 'LineWidth', 2); hold on;
title('Cosine Gabor Filter');
xlabel('Time');
ylabel('Amplitude');
grid on;

subplot(1, 3, 2);
plot(x, sine_gabor_shifted, 'LineWidth', 2);
title('Sine Gabor Filter');
xlabel('Time');
ylabel('Amplitude');
grid on;

subplot(1, 3, 3);
plot(x, sine_sine_gabor, 'LineWidth', 2);
title('Derivaive of Sine Gabor Filter');
xlabel('Time');
ylabel('Amplitude');
grid on;

% Function to compute Gabor filter
function gb = gabor_filter(x, frequency, sigma, phase)
    if nargin < 4
        phase = 0; % Default phase is 0 for cosine Gabor filter
    end
    gb = exp(-0.5 * (x .^ 2 / sigma ^ 2)) .* cos(2 * pi * frequency * x + phase);
end
