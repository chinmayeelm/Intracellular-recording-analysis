
function [raster_data,avg_gcfr,no_of_true_trials, gcfr, invalid_trials]   = get_raster_gcfr(no_of_trials, P_rec, single_trial_length, fs, L, alpha)
    
        gcfr = [];
        invalid_trials = [];
        no_of_true_trials = 0;
        raster_data = zeros(no_of_trials, single_trial_length);
%         L = fs/2; %fs/10; %1000; %5000
%         alpha = 4;%2;%8;
        gauss_win = gausswin(L, alpha); 
        for i=1:no_of_trials
%             p=[]; l=[];
            [p,l] =  findpeaks(P_rec(i,:), "MinPeakHeight",0.3*max(P_rec(i,:)), "MinPeakDistance", 0.003*fs);
            spike_amp= p;
            ISI = diff(l)./fs;
%             plot(spike_amp(2:end),ISI, '.'); hold on;
%             figure;
%             scatter(spike_amp(2:end),ISI, '.'); %hold on;
%             plot(i, spike_amp, '.'); hold on;
%             ylabel 'Amplitude (mV)'
%             xlabel 'Trial No.'
            
%             plot(i, ISI, '.'); hold on;
%             ylabel 'ISI (ms)';
%             xlabel 'Trial No.';
%             ylim([0 1]);
%             xlim([0 47]);
            if mode(p)<5 %3 
                disp("trial=");disp(i);
                disp("mode(p)="); disp(mode(p));
                disp("min(rec)="); disp(min(P_rec(i,:)));
%                 disp("p="); disp(p);
                raster_data(i,:) = nan;
                invalid_trials = [invalid_trials; i];
                continue;
            else
                raster_data(i,l) = 1;
                no_of_true_trials = no_of_true_trials+1;
                gcfr(no_of_true_trials,:) = (fs/L)*filter(gauss_win, 1, raster_data(i,:));
%                 gcfr(no_of_true_trials,:) = filtfilt(gauss_win, 1, raster_data(i,:));

            end
        end
        
        raster_data(invalid_trials,:)=[];
        sum_of_spikes = sum(raster_data, 1);
        
        avg_gcfr = (fs/L)*(filter(gauss_win, 1, sum_of_spikes))/no_of_true_trials;

end    