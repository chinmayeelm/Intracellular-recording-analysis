function [timeBinCenters,firingRate] = psth(raster,binWidth, fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    timeBins = 1 : binWidth*fs : length(raster)+1;
    for iBin =  1:length(timeBins)-1
        firingRate(iBin) = sum(raster(:,timeBins(iBin):timeBins(iBin+1)-1),"all")/(binWidth*size(raster,1));
        timeBinCenters(iBin) = (timeBins(iBin)+(binWidth*fs/2))/fs;
   
    end
end