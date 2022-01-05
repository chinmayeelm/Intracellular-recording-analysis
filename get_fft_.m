function [stim_freq, power_fft, frq_fft] = get_fft_(stim_matrix, Fs, L)

Y = fft(stim_matrix)/L;

NyLimit = Fs/2;

F = linspace(0,1,L/2)*NyLimit;

figure;
plot(F,abs(Y(1:L/2)))

power_fft = Y;
frq_fft = F;

[~, loc] = max(power_fft(2:end));
stim_freq = F(loc+1);


end