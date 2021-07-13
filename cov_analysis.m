function [cov_matrix, ev1, ev2] = cov_analysis(raster_data, stimulus, window, fs)

    [m,~] = size(raster_data);
    all_spike_triggers = [];
    for i=1
        spike_locs = find(raster_data(i,:)==1);
        if isempty(spike_locs)
            continue;
        end
        for j=1:100 %length(spike_locs)
            if (spike_locs(j)-window*fs)<= 0 
                continue;
            end
            spike_triggers(j,:) = stimulus(i,(spike_locs(j)-window*fs):spike_locs(j));
            spike_triggers(j,:) = spike_triggers(j,:)-mean(spike_triggers(j,:));
 
        end

        all_spike_triggers = [all_spike_triggers; spike_triggers];

    end
    pattern_length = window*fs+1;
    stimulus_prior = reshape(stimulus(1,1:100*pattern_length), [100,pattern_length]);
    size(stimulus_prior)
    size(all_spike_triggers)
    nsp = size(all_spike_triggers,1); %number of spikes or no. of stimulus patterns
    STA = mean(all_spike_triggers,1);
    size(STA)
    
%     cov_matrix = (1/(length(all_spike_triggers)-1)).*((all_spike_triggers-STA)'*(all_spike_triggers-STA));
%     cov_matrix = (1/(length(all_spike_triggers)-1)).*(all_spike_triggers*all_spike_triggers');
    
    cov_matrix = cov((all_spike_triggers'-STA'));
    figure;
    heatmap(flip(cov_matrix,2), 'Colormap', jet);
    title('Covariance matrix using the formula');
    ylabel('spike triggering stimulus pattern');
    xlabel('spike triggering stimulus pattern');
%     figure;
%     heatmap(flip(cov_matrix_direct,2), 'Colormap', jet);
%     title('Covariance matrix using "cov" command');
%     ylabel('spike triggering stimulus pattern');
%     xlabel('spike triggering stimulus pattern');
    
%     figure;
%     heatmap(flip(cov_matrix-cov_matrix_direct,2));

%     size(cov_matrix)

%     stimulus_covariance = (1/(length(stimulus_prior)-1)).*((stimulus_prior)'*(stimulus_prior));
    stimulus_covariance = cov(stimulus_prior');
    figure;
    heatmap(flip(stimulus_covariance,2), 'Colormap', jet);
    title('Stimulus covariance using the formula');
    ylabel('stimulus pattern');
    xlabel('stimulus pattern');
    
    diff_cov = cov_matrix - stimulus_covariance;
    figure;
    heatmap(flip(cov_matrix,2), 'Colormap', jet);
    title('STC-C');
%     ylabel('spike triggering stimulus pattern');
%     xlabel('spike triggering stimulus pattern');
    
    [V,D] = eig(diff_cov);
    
    [d, ind] = sort(diag(D));
    
    D_sorted = D(ind, ind);
    
    V_sorted = V(:,ind);
    
    ev1 = V_sorted(1,:);
    ev2 = V_sorted(2,:);
    ev_last = V_sorted(end,:);
    size(V_sorted)
    
    figure;
    t = linspace(-window*1000, 0,length(ev1));
    plot(t, V_sorted(1:2,:));
%     plot(t, ev1, t, ev2, t, ev_last);
%     plot(ev1);
%     figure; plot(t,ev2);
end




