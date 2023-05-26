cd 'E:\Recordings\LUTs';
% rampData = readlines('validRampList.txt');
% ignorefiles = [1 2 4 5 6 8 17 19 25 27 30 33 42 44 46];
% rampData(ignorefiles, :) = [];
% rampData = readlines('rampDataList.txt');
% dataDirectory_ramp = extractBefore(rampData, '_');
% filename_ramp = extractAfter(rampData, "_");
% nfiles = length(filename_ramp);

stepData = readlines('allStepList.txt');
dataDirectory_step = extractBefore(stepData, '_');
filename_step = extractAfter(stepData, "_");
nfiles = length(filename_step);

% stairData = readlines('allStairList.txt');
% dataDirectory_stair = extractBefore(stairData, '_');
% filename_stair = extractAfter(stairData, "_");
% nfiles = length(filename_stair);

% adaptation_coeff = zeros(1, nfiles);
% rsquare = zeros(1, nfiles);

ignorefiles = [];

for irow = 21%:nfiles
    
    irow
    close all;
    
%     try
%     [P_ramp,intendedStimulus] = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),0);
    
    [P_step, intendedStimulus] = getStructP(dataDirectory_step(irow), filename_step(irow),0);
%     [P_stair, intendedStimulus] = getStructP(dataDirectory_stair(irow), filename_stair(irow),0);
%     rampGCFRplots(P_ramp, intendedStimulus);
    stepGCFRplots(P_step);
%     stairGCFRplots(P_stair, intendedStimulus);
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