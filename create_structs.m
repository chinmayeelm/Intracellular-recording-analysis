function P = create_structs(rec_protocols_sorted,stim_protocols_hes_sorted,fs, stim_protocols_ifb_sorted, no_of_protocols, no_of_trials, single_trial_length, stim_order_sorted, max_chirp_frq,amp_sweep_frq, blwgn_fc, max_sqr_sin_frq, impulse_dur)

    for i=1:no_of_protocols
        
        
        P(i).rec = 100*( rec_protocols_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:)); %voltages will be in mV
        P(i).hes_data_unfilt = stim_protocols_hes_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        P(i).stim_ifb = stim_protocols_ifb_sorted(1+(i-1)*no_of_trials:no_of_trials*i,:);
        
        P(i).stim_name = stim_order_sorted(i*no_of_trials);
        type_frq = split(P(i).stim_name, '_');
        P(i).stim_type = type_frq(1);
        
        
        if (P(i).stim_type == "sin" ) || (P(i).stim_type == "sqr")
            P(i).stim_period = 1/str2num(type_frq(2));
            order=4;
            for j=1:no_of_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 1.5*str2double(type_frq(2)), fs, order);
            end
            
        elseif P(i).stim_type == "frq" || P(i).stim_type == "dec"
            order=8;
            for j=1:no_of_trials
                P(i).antennal_movement(j, :) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), max_chirp_frq, fs, order);
                P(i).max_chirp_frq = max_chirp_frq;
            end
            
        elseif P(i).stim_type == "impulse"
            P(i).imp_dur = impulse_dur;
            order=4;
            for j=1:no_of_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 1/P(i).imp_dur, fs, order);
            end
            
        elseif P(i).stim_type == "blwgn"
            order=4;
            for j=1:no_of_trials
                P(i).antennal_movement(j,:) = butter_filtfilt(P(i).hes_data_unfilt(j,:), blwgn_fc, fs, order);
                P(i).blwgn_fc = blwgn_fc;
            end
            
        elseif (P(i).stim_type == "amp")
            order=3;
            for j=1:no_of_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), amp_sweep_frq, fs, order);
                P(i).amp_sweep_frq = amp_sweep_frq; 
                P(i).stim_period = 1/amp_sweep_frq;
            end
        
        elseif (P(i).stim_type == "step")
            P(i).step_dur = type_frq(2);
            order=4;
            for j=1:no_of_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2/str2double(type_frq(2)), fs, order);
            end
            
        elseif (P(i).stim_type == "sum")
            P(i).stim_period = 1/str2num(type_frq(3));
            order=4;
            for j=1:no_of_trials
                P(i).antennal_movement(j,:) =  butter_filtfilt(P(i).hes_data_unfilt(j,:), 2*str2double(type_frq(3)), fs, order);
            end
            
            
        end

        
        
        [raster,avg_gcfr,complete_trials, gcfr]  = get_raster_gcfr(no_of_trials, P(i).rec, single_trial_length);
        P(i).gcfr = gcfr;
        P(i).raster = raster;
        P(i).avg_gcfr = avg_gcfr;
        P(i).complete_trials = complete_trials;

    end

end