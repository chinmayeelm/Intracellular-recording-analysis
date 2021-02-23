function [stim_freq, power_fft, frq_fft] = fft_stim(stim_matrix, Fs, L)


%     [m,~]  = size(stim_matrix);
%     stim_freq = zeros(m, 1);
    
%     for i=1:m
%         stim_filt = sgolayfilt(stim_matrix(i,:), 3, 51);
        Y = fft(stim_matrix);

%         figure()
%         plot(stim_matrix(i, :));

        P2 = abs(Y/L);
        P1 = P2(1:(L/2)+1);
        P1(2:end-1) = P1(2:end-1);
%         length(P1(2:end-1))

       
        f = Fs*(0:(L/2))/L;
%         f = f(2:end-1);
        length(f)
        figure();
%         plot(f,P1(2:end-1)); hold on;
        plot(f,P1); %hold on;
        title('Tuning curve')
        xlabel('f (Hz)')
        ylabel('|P1(f)|')
        xlim([0 350]);

%         [~,loc] = max(P1(2:end-1));
        [~,loc] = max(P1);
        stim_freq = f(loc);
        
        power_fft = P1;
        frq_fft = f;
    
%     end

end