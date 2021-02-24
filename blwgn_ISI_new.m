ISI_all = [];
for i =1
    
    fs = blwgn(i).fs;
    start = blwgn(i).OFF_dur *fs;
    stop = (blwgn(i).OFF_dur + blwgn(i).ON_dur) *fs;
    raster = blwgn(i).raster(:,start:stop);
    [rows,~] = size(raster);
    factor = 1000/fs; %1000 ms/fs to convert ISI to ms
    
    for j=1:rows
        locs = find(raster(j,:)==1);
        
        
        ISI = diff(locs,1,2).*factor;
        t = linspace(blwgn(i).OFF_dur,(blwgn(i).OFF_dur + blwgn(i).ON_dur), length(ISI)) ;
    
    
        plot(ISI, '.'), hold on;
        
        ISI_all = [ISI_all ISI];

%         figure;
%         scatter(t, ISI); hold on;
    end
     
    ylabel 'ISI (ms)';
    xlabel 'time (s)'  
    
end

figure;
histogram(ISI_all,'BinWidth', 0.001);
xlabel 'ISI (ms)';
ylabel 'Occurances';
title 'Histogram of ISI';