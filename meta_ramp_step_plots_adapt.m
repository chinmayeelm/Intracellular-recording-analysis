cd 'E:\Recordings\LUTs';
rampData = readlines('validRampList.txt');
% ignorefiles = [1 2 4 5 6 8 17 19 25 27 30 33 42 44 46];
% rampData(ignorefiles, :) = [];
% rampData = readlines('rampDataList.txt');
dataDirectory_ramp = extractBefore(rampData, '_');
filename_ramp = extractAfter(rampData, "_");
nfiles = length(filename_ramp);

% stepData = readlines('allStepList.txt');
% dataDirectory_step = extractBefore(stepData, '_');
% filename_step = extractAfter(stepData, "_");
% nfiles = length(filename_step);

% stairData = readlines('allStairList.txt');
% dataDirectory_stair = extractBefore(stairData, '_');
% filename_stair = extractAfter(stairData, "_");
% nfiles = length(filename_stair);

% chirpData = readlines('chirpPartialList.txt');
% dataDirectory_chirp = extractBefore(chirpData, '_');
% filename_chirp = extractAfter(chirpData, "_");
% nfiles = length(filename_chirp);

% adaptation_coeff = zeros(1, nfiles);
% rsquare = zeros(1, nfiles);

ignorefiles = [];

for irow = 1:nfiles
    
    irow
    close all;
    
%     try
    P_ramp = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),0);
%     P_step = getStructP(dataDirectory_step(irow), filename_step(irow),0);
%     P_stair = getStructP(dataDirectory_stair(irow), filename_stair(irow),0);
%     P_chirp = getStructP(dataDirectory_chirp(irow), filename_chirp(irow),0);
    
    rampGCFRplots(P_ramp);
%     stepGCFRplots(P_step);
%     stairGCFRplots(P_stair);
%     chirpGCFRplots(P_chirp);
    try
    [adaptation_coeff(irow), rsquare(irow)] = calcAdaptation(P_ramp);
    catch
        disp("Not enough duration to fit curve");
        ignorefiles = [ignorefiles irow];
    end
    
%     catch
%         ignorefiles = [ignorefiles irow];
%     end
    pause;

end