fs = M1N1_chirp.fs;

cmp = summer(10);
num = [1 2]; %[3 17:20 23:27];
% str = strings(1,length(num));
% str(:) = "T";
% struct_list = strcat(str',num2str(num'));

%%
for i=1:length(num)
%     P = M1N1_chirp.(struct_list(i)); 
    
    ON_dur = M1N1_chirp.ON_dur;
    OFF_dur = M1N1_chirp.OFF_dur;

    % zero_crossing(M1N1_chirp.stim_ifb(1,:), M1N1_chirp.norm_gcfr, fs, ON_dur, OFF_dur);
    [stim_freq, max_FR] = tuning_curve(M1N1_chirp(i).antennal_movement, M1N1_chirp(i).norm_gcfr, fs, ON_dur, OFF_dur); %require edit in the function to take all rows
    figure(1+mod(i,5));
%     figure();
    [lineOut, fillOut] = stdshade(max_FR,0.3,cmp(i,:),stim_freq); hold on;
    lineOut.LineWidth = 0.3;
    ylabel('Normalised firing rate (spike/s)');
    xlabel('Stimulus frequency (Hz)');
%     title(strcat("Chirp duration = ", num2str(M1N1_chirp.ON_dur), " s"));
    xlim([0 120]);
%     lgd = legend('','Increasing chirp','', 'Decreasing chirp');
%     lgd.Location = 'northwest';


    [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(M1N1_chirp(i).antennal_movement(1,:), M1N1_chirp(i).raster(1,:), fs, ON_dur, OFF_dur);
    figure(6+mod(i,5)); scatter(I_spike_freq,I_spike_phase,100, 'k.'); hold on; 
    pmin = min(I_spike_phase);
    pmax = max(I_spike_phase);
    pimin = floor(pmin/pi);
    pimax = ceil(pmax/pi);
    %              yticks(0:pi/4:2*pi);
    yticks((pimin:pimax) * pi);
    yticklabels( string(pimin:pimax) + "\pi" )
    scatter(II_spike_freq,II_spike_phase, 100,  'b.');  
    scatter(III_spike_freq,III_spike_phase, 100, 'r.')
    ylabel('Spike phase (rad)');
    xlabel('Stimulus frequency (Hz)');
    legend('I spike','II spike', 'III spike');
    title(strcat("Chirp duration = ", num2str(M1N1_chirp.ON_dur), " s"));

%     M1N1_chirp.frq_chirp_f = stim_freq;
%     M1N1_chirp.frq_chirp_FR = max_FR;
% 
%     M1N1_chirp.I_spike = [I_spike_freq', I_spike_phase'];
%     M1N1_chirp.II_spike = [II_spike_freq', II_spike_phase'];
%     M1N1_chirp.III_spike = [III_spike_freq', III_spike_phase'];

end