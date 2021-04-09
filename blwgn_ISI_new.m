
for i =1:9
    ISI_all = [];
    fs = blwgn(i).fs;
    start = blwgn(i).OFF_dur *fs;
    stop = (blwgn(i).OFF_dur + blwgn(i).ON_dur) *fs;
    raster = blwgn(i).raster(:,start:stop);
%     rec = blwgn(i).rec(:,start:stop);
%     stim = blwgn(i).antennal_movement(:,start:stop);
%     gcfr = blwgn(i).norm_gcfr(:,start:stop);
    [rows,~] = size(raster);
    factor = 1000/fs; %1000 ms/fs to convert ISI to ms
%     spike_amp = struct();
%     isi = struct();
%     delta_amp_all = [];
    for j=1:rows
        locs = find(raster(j,:)==1);
        
%         [spike_amp,locs] =  findpeaks(rec(j,:), "MinPeakHeight",0.2*max(rec(j,:)));
        
        ISI = diff(locs,1,2).*factor;  %ISI in ms
%         delta_amp = diff(spike_amp, 1, 2);
        t = linspace(0, blwgn(i).ON_dur, length(raster)) ;
        
%         Vm = rec(j,:);
%         diff_Vm = diff(Vm)/0.1;
        
%         diff_stim = diff(stim(j,:))/0.0001;
%         figure;
%         A1 = subplot(2,1,2); plot(t(2:end), diff_stim); hold on;
%         figure;
%         A2 = subplot(2,1,1); plot(locs(2:end)/fs, ISI, 'r.');
%         A2 = subplot(2,1,1); plot(t, Vm);
        
%         linkaxes([A1 A2], 'x');
        
%         figure;
%         plot(t, rec(j,:), locs/fs, spike_amp, 'r.');
%         figure;
%         plot(ISI, '.'), hold on;
%         plot(ISI, spike_amp(2:end), '.'); hold on;

%             figure;
%         scatter(ISI, delta_amp, '.'); hold on;
%         plot(j, ISI, '.'); hold on;
%         subplot(20,1,j); plot(t, gcfr(j,:), 'k');
%         xlabel 'ISI (ms)'
%         ylabel 'Spike amplitude'
%         ylim ([-10 0]);
%         axis off
%         xlim ([0 21])
%             ylim([0 1]);
%             xlim([0 47]);
        
        ISI_all = [ISI_all ISI];
%         delta_amp_all = [delta_amp_all delta_amp];

%         figure;
%         scatter(t, ISI); hold on;
    end
     
%     ylabel 'ISI (ms)';
%     xlabel 'time (s)'  
%     figure;
    subplot(3,3,i);
    histogram(ISI_all, 'BinWidth', .001);
    xlabel 'ISI (ms)';
    ylabel 'Occurances';
    title (strcat('Histogram of ISI',"   ", blwgn(i).filename), 'Interpreter', 'none');
    xlim ([5 45]);
    ylim ([0 250]);
    
    filename = strcat("histISI_", blwgn(i).filename);
    saveas(gcf,filename,'png');
end

% l = find(delta_amp_all<0);
% neg_delta_amp = delta_amp_all(l);
% neg_del_amp_isi = ISI_all(l);
% figure;
% scatter(neg_del_amp_isi,neg_delta_amp, '.');
