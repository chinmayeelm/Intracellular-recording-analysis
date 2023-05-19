rampData = readlines('allRampList.txt');
% ignorefiles = [1 2 4 5 6 8 17 19 25 27 30 33 42 44 46];
% rampData(ignorefiles, :) = [];
% rampData = readlines('rampDataList.txt');
stepData = readlines('stepDataList.txt');
dataDirectory_ramp = extractBefore(rampData, '_');
filename_ramp = extractAfter(rampData, "_");
dataDirectory_step = extractBefore(stepData, '_');
filename_step = extractAfter(stepData, "_");
nfiles = length(filename_ramp);
adaptation_coeff = zeros(1, nfiles);
rsquare = zeros(1, nfiles);

ignorefiles = [];

for irow = 6%:nfiles
    
    irow
%     close all;
    
%     try
    P_ramp = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),0);
    
%     P_step = getStructP(dataDirectory_step(irow), filename_step(irow),0);
    rampGCFRplots(P_ramp);
%     stepGCFRplots(P_step);
%     try
%     [adaptation_coeff(irow), rsquare(irow)] = calcAdaptation(P_ramp);
%     catch
%         disp("Not enough duration to fit curve");
%         ignorefiles = [ignorefiles irow];
%     end
    
%     catch
%         ignorefiles = [ignorefiles irow];
%     end
    pause;

end