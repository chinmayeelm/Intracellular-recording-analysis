ON_dur = 10;
OFF_dur = 5; 
fs = 10000;

start = fs*OFF_dur; stop = (ON_dur+OFF_dur)*fs;

stim = P.antennal_movement(:, start:stop);
resp = P.raster(:, start:stop);
gcfr = P.gcfr(:,start:stop);
[rows,~] = size(resp);
isi_all = [];

    for i=1:rows
        spike_locs = find(resp(i,:));
        isi = (diff(spike_locs)*1000)/fs;
        isi_all = [isi_all isi];
    end
    
figure;
histogram(isi,'BinWidth', 0.001);
xlabel 'ISI (ms)';
ylabel 'Occurances';
title 'Histogram of ISI';
