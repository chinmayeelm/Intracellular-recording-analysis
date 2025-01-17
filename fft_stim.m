function [stim_freq, power_fft, frq_fft] = fft_stim(stim_matrix, Fs)


%     [m,~]  = size(stim_matrix);
%     stim_freq = zeros(m, 1);
    
%     for i=1:m
%         stim_filt = sgolayfilt(stim_matrix(i,:), 3, 51);
        L = length(stim_matrix);
        n = L;%2^nextpow2(L);
        
        Y = fft(stim_matrix,n);

%         figure()
        % plot(stim_matrix);

        P2 = abs(Y/n);
        P1 = P2(1:floor(n/2+1));
        P1(2:end-1) = 2*P1(2:end-1);
        length(P1(2:end-1));

       
        % f = Fs*(0:(L/2))/L;
        f = Fs*(0:(n/2))/n;
%         f = f(2:end-1);
        length(f);
%         figure();
        % yyaxis right; 
        % plot(f(2:end-1),P1(2:end-1), 'LineWidth',1); hold on;
        % % plot(f,P1); %hold on;
        % title('FFT of Stimulus')
        % xlabel('f (Hz)')
        % ylabel('Power')
%         xlim([0 350]);

%         [~,loc] = max(P1(2:end-1));
%         [~,loc] = max(P1);



        [~,loc] = max(P1);
        stim_freq = f(loc);
        
        power_fft = P1(2:end-1);
        frq_fft = f(2:end-1);
    
%     end

end