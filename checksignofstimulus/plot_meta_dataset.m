latency = [];
k=1;
for i=1:length(meta_struct)
    for j=1:length(meta_struct(i).moth)
        if (meta_struct(i).moth(j).stim_type=="sqr") && (meta_struct(i).moth(j).stim_period == 8)
            
%               h = figure(); %subplot(3, 1,k);
%                 
% %               FR_normalised = meta_struct(i).moth(j).frq_chirp_FR./max(meta_struct(i).moth(j).frq_chirp_FR);
%               [lineOut, ~] = stdshade(meta_struct(i).moth(j).frq_chirp_FR, 0.5, [0.4940 0.1840 0.5560], meta_struct(i).moth(j).frq_chirp_f);
%               %h.XAxis.Visible = 'off';
% %               h.XLim = [50 300];
% %               h.YLim = [0 1.2];
%               box off;
%               k=k+1;

%               fig = figure();
%               t = meta_struct(i).moth(j).time(1:meta_struct(i).moth(j).single_trial_length);
%               subplot(2,1,1); stdshade(meta_struct(i).moth(j).norm_gcfr,0.5, [0.4660 0.6740 0.1880],t);
%               ylabel('Normalised GCFR');
%               subplot(2,1,2); stdshade(meta_struct(i).moth(j).antennal_movement,0.5, [0.6350 0.0780 0.1840],t);
%               ylabel('Antennal movement (mm)');
%               xlabel('time (s)');
%             title((join(split(meta_struct(i).moth(j).stim_name,"_")," ")) +" Hz");
%             
%             savefigures(meta_struct(i).moth(j).filename, meta_struct(i).moth(j).stim_name, fig, meta_struct(i).moth(j).date);
%                 t= -100:1/10:0;
%               plot(t, meta_struct(i).moth(j).STA), hold on;
%               ylabel('Antennal movement (mm)');
%               xlabel('time before spike (ms)');
%               box off;
%               colormap(f, jet);

%                 plot(meta_struct(i).moth(j).frq_fft, meta_struct(i).moth(j).power_fft./max(meta_struct(i).moth(j).power_fft)); hold on;
%                 ylabel('Normalised STA response');
%                 xlabel('Frequency (Hz)');
%                 xlim([0 350]);
%                 box off;
%                 
                
              
                

        end
    end
end
% colormap(gcf, parula);
% figure();
% histogram(latency, 250);
% length(latency)
% ylabel('Occurences'); xlabel('Latency (s)');

