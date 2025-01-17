function P = create_structs(rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted,intended_stimulus_sorted, no_of_protocols, no_of_trials, single_trial_length, stim_order_sorted, max_chirp_frq,amp_sweep_frq, blwgn_fc, impulse_dur, a, b, c, L, sigma, movementRadius)

    for i=1:no_of_protocols
        
        
        
        stim_names = unique(stim_order_sorted);
%         P(i).stim_name = stim_order_sorted(i*no_of_trials);
        P(i).stim_name = stim_names(i);
        
        
        if length(stim_order_sorted) == 1 && no_of_trials > 1
            no_of_trials;
        else
            no_of_trials = length(find(stim_order_sorted == P(i).stim_name));
        end
        
%         P(i).rec = 100*( rec_protocols_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:)); %voltages will be in mV
%         P(i).hes_data_unfilt = stim_protocols_hes_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
%         P(i).stim_ifb = stim_protocols_ifb_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        
        trial_idx = find(stim_order_sorted == P(i).stim_name,1);
        [trial_idx  trial_idx+ no_of_trials-1];
        P(i).rec = 100*(rec_protocols_sorted(trial_idx: trial_idx+ no_of_trials-1, :));
        P(i).hes_data_unfilt = stim_protocols_hes_sorted(trial_idx: trial_idx+ no_of_trials-1, :);
        P(i).stim_ifb = stim_protocols_ifb_sorted(trial_idx: trial_idx+ no_of_trials-1, :);
        P(i).intendedStimulus = intended_stimulus_sorted(trial_idx: trial_idx+ no_of_trials-1, :);
        
        
        [raster,avg_gcfr,complete_trials, gcfr, invalid_trials]  = get_raster_gcfr(no_of_trials, P(i).rec, single_trial_length, fs, L, sigma);
        P(i).gcfr = gcfr; %in Hz
        P(i).raster = raster;
        P(i).avg_gcfr = avg_gcfr; %in Hz
        P(i).complete_trials = complete_trials;
        
        P(i).hes_data_unfilt(invalid_trials,:) = [];
        P(i).stim_ifb(invalid_trials, :) = [];
        P(i).rec(invalid_trials,:) = [];
        P(i).intendedStimulus(invalid_trials,:) = [];
        
        type_frq = split(P(i).stim_name, '_');
        P(i).stim_type = type_frq(1);
        
        
        if (P(i).stim_type == "sin" )
            P(i).stim_period = 1/str2double(type_frq(2));
            order=4;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 10*str2double(type_frq(2)), fs, order, a, b, c); %2*str2double(type_frq(2))
%              
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif (P(i).stim_type == "sqr")
            P(i).stim_period = 1/str2double(type_frq(2));
            order=4;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 10, fs, order, a, b, c); %2*str2double(type_frq(2))
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif P(i).stim_type == "frq" || P(i).stim_type == "dec"
            order=4;
            for j=1:complete_trials
                P(i).antennal_movement(j, :) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2*max_chirp_frq, fs, order, a, b, c);
                P(i).max_chirp_frq = max_chirp_frq;
            end
            % P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - P(i).antennal_movement(1,1))./movementRadius);
            
        elseif P(i).stim_type == "impulse"
            P(i).imp_dur = impulse_dur;
            order=4;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2/P(i).imp_dur, fs, order, a, b, c);
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif P(i).stim_type == "blwgn" || P(i).stim_type == "blwgn2"
            order=10;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) = butter_filtfilt(P(i).hes_data_unfilt(j,:), blwgn_fc, fs, order, a, b, c);
                P(i).blwgn_fc = blwgn_fc;
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif (P(i).stim_type == "amp")
            order=3;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2*amp_sweep_frq, fs, order, a, b, c);
                P(i).amp_sweep_frq = amp_sweep_frq; 
                P(i).stim_period = 1/amp_sweep_frq;
            end
            % P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif (P(i).stim_type == "step" || type_frq(2) == "ramp" || P(i).stim_type == "var" )
%             P(i).step_dur = type_frq(2);
            order=4;
            for j=1:complete_trials
%                 P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2/str2double(type_frq(2)), fs, order, a, b, c);
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 20, fs, order, a, b, c); %20 % 8 before 30.11.2023
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);

        elseif ( type_frq(2) == "stair" || P(i).stim_name == "inc_dec_stair")
%             P(i).step_dur = type_frq(2);
            order=4;
            for j=1:complete_trials
%                 P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2/str2double(type_frq(2)), fs, order, a, b, c);
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 20, fs, order, a, b, c); % fc was 10 earlier
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif (P(i).stim_type == "sum")
            P(i).stim_period = 1/str2num(type_frq(3));
            order=4;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 4*str2double(type_frq(3)), fs, order, a, b, c);
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        elseif (P(i).stim_type == "noisySin")
            noisySin_frq = 30;
            max_noise_frq = 300;
            order = 4;
            for j=1:complete_trials
                P(i).antennal_movement(j,:) = butter_filtfilt(P(i).hes_data_unfilt(j,:), 2*max_noise_frq, fs, order, a, b, c);
            end
            P(i).antennal_movement = rad2deg((P(i).antennal_movement - mean(P(i).antennal_movement(:,1:2*fs),2))./movementRadius);
            
        end

        mean_movement = mean(P(i).antennal_movement,1);
%         disp("protocol="); disp(P(i).stim_name);
%         disp("Max. antennal movement in mm ="); disp(max(mean_movement)-min(mean_movement));
        
        P(i).mean_movement = mean_movement;

    end

end