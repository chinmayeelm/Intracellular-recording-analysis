function [raster_data,avg_gcfr,no_of_true_trials, gcfr]   = get_raster_gcfr(no_of_trials, P_rec, single_trial_length)
    

        no_of_true_trials = 0;
        raster_data = zeros(no_of_trials, single_trial_length);
        L = 1000;
        alpha = 4;
        gauss_win = gausswin(L, alpha);
        for i=1:no_of_trials
%             p=[]; l=[];
            [p,l] =  findpeaks(P_rec(i,:), "MinPeakHeight",0.2*max(P_rec(i,:)));

            if mode(p)<2%5
                continue;
            else
                raster_data(i,l) = 1;
                gcfr(i,:) = filter(gauss_win, 1, raster_data(i,:));
                no_of_true_trials = no_of_true_trials+1;
            end
        end   

        sum_of_spikes = sum(raster_data, 1);

        
        avg_gcfr = (filter(gauss_win, 1, sum_of_spikes))/no_of_true_trials;

end    