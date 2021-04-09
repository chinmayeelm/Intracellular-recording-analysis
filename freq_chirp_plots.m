fs = frq_chirp(1).fs;

cmp = summer(10);
num = [1 2]; %[3 17:20 23:27];
% str = strings(1,length(num));
% str(:) = "T";
% struct_list = strcat(str',num2str(num'));

%%
for i=1
%     P = frq_chirp.(struct_list(i)); 
    
    ON_dur = frq_chirp(i).ON_dur;
    OFF_dur = frq_chirp(i).OFF_dur;

%     figure(1+mod(i,5));
    figure();
    start = OFF_dur*fs;
    stop = (ON_dur+OFF_dur)*fs;
%     stim_freq = linspace(1,frq_chirp(i).max_chirp_frq,ON_dur*fs+1); %inc chirp
    stim_freq = linspace(1,frq_chirp(i).max_chirp_frq,ON_dur*fs+1); %dec chirp
    [lineOut, fillOut] = stdshade(frq_chirp(i).norm_gcfr(:,start:stop),0.3,cmp(i,:),stim_freq); hold on;
    lineOut.LineWidth = 0.3;
    ylabel('Normalised firing rate (spike/s)');
    xlabel('Stimulus frequency (Hz)');
    title(strcat("Chirp duration = ", num2str(ON_dur), " s"));
%     xlim([0 300]);
%     lgd = legend('','Increasing chirp','', 'Decreasing chirp');
%     lgd.Location = 'northwest';


%     [I_spike_phase, II_spike_phase, III_spike_phase, I_spike_freq, II_spike_freq, III_spike_freq] = spike_phase(frq_chirp(i).antennal_movement(1,:), frq_chirp(i).raster(1,:), fs, ON_dur, OFF_dur);
% %     figure(6+mod(i,5)); 
%     figure;
%     scatter(I_spike_freq,I_spike_phase,100, 'k.'); hold on; 
%     pmin = min(I_spike_phase);
%     pmax = max(I_spike_phase);
%     pimin = floor(pmin/pi);
%     pimax = ceil(pmax/pi);
%     %              yticks(0:pi/4:2*pi);
%     yticks((pimin:pimax) * pi);
%     yticklabels( string(pimin:pimax) + "\pi" )
%     scatter(II_spike_freq,II_spike_phase, 100,  'b.');  
%     scatter(III_spike_freq,III_spike_phase, 100, 'r.')
%     ylabel('Spike phase (rad)');
%     xlabel('Stimulus frequency (Hz)');
%     legend('I spike','II spike', 'III spike');
%     title(strcat("Chirp duration = ", num2str(ON_dur), " s"));

%     frq_chirp.dec_frq_chirp_f = stim_freq;
%     frq_chirp.dec_frq_chirp_FR = max_FR;
% 
%     frq_chirp.I_spike = [I_spike_freq', I_spike_phase'];
%     frq_chirp.II_spike = [II_spike_freq', II_spike_phase'];
%     frq_chirp.III_spike = [III_spike_freq', III_spike_phase'];

end