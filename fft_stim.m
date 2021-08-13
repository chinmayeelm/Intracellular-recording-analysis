function [power_fft, frq_fft] = fft_stim(stim_matrix, Fs, N)

    L = length(stim_matrix);
    dF = Fs/L;                      % hertz
    f = -Fs/2:dF:Fs/2-dF;
    stim_fft = fftshift(fft(stim_matrix));

%     figure;
%     plot(f,abs(stim_fft)/L); hold on;
% 
%     xlabel('Frequency (in hertz)');
%     ylabel('Response magnitude');
%     xlim ([0 350]);
    
   
    power_fft = abs(stim_fft(length(stim_fft)/2 +1 :end));
    frq_fft = f(length(f)/2 +1:end);

end