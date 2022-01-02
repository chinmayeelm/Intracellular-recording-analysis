function isi_all = isi_hist(resp, fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[rows,~] = size(resp);
isi_all = [];

    for i=1:rows
        spike_locs = find(resp(i,:));
        isi = (diff(spike_locs)*1000)/fs;
        isi_all = [isi_all isi];
    end
    
% figure;
% histogram(isi_all,'BinWidth', 1);
% xlabel 'ISI (ms)';
% ylabel 'Occurances';
% title 'Histogram of ISI';

end

