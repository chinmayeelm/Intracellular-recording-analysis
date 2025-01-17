for i=1:9
plot(blwgn(i).frq_fft, blwgn(i).power_fft/max(blwgn(i).power_fft)); hold on;
xlim ([0 350]);
end

ylabel 'Power';
xlabel 'Frequency (Hz)';
title ("tuning curves");
