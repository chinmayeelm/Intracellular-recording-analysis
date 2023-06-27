cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
% rampData = readlines('validRampList.txt');
% dataDirectory_ramp = extractBefore(rampData, '_');
% filename_ramp = extractAfter(rampData, "_");
% nfiles = length(filename_ramp);

% stepData = readlines('stepDataList.txt');
% dataDirectory_step = extractBefore(stepData, '_');
% filename_step = extractAfter(stepData, "_");
% nfiles = length(filename_step);

stairData = readlines('allStairList.txt');
dataDirectory_stair = extractBefore(stairData, '_');
filename_stair = extractAfter(stairData, "_");
nfiles = length(filename_stair);

% chirpData = readlines('chirpPartialList.txt');
% dataDirectory_chirp = extractBefore(chirpData, '_');
% filename_chirp = extractAfter(chirpData, "_");
% nfiles = length(filename_chirp);

% wnData = readlines('blwgnValidData.txt');
% dataDirectory_wn = extractBefore(wnData, '_');
% filename_wn = extractAfter(wnData, "_");
% nfiles = length(filename_wn);

% adaptation_coeff = zeros(1, nfiles);
% rsquare = zeros(1, nfiles);

ignorefiles = [];
velocity = [];
FR = [];

% figure;
c = parula(21);
% ci = 0;
for irow = 15%1:nfiles
    
    irow
%     close all;
    ci = ci+1;
%     try
    % P_ramp = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),0);
    % P_step = getStructP(dataDirectory_step(irow), filename_step(irow),0);
    P_stair = getStructP(dataDirectory_stair(irow), filename_stair(irow),0);
    % P_chirp = getStructP(dataDirectory_chirp(irow), filename_chirp(irow),0);
%     P_wn = getStructP(dataDirectory_wn(irow), filename_wn(irow),0);
    % title(irow)

%     S(irow) = jitter_func(P_wn, irow);
    % [vel, meanMaxFR] = rampGCFRplots(P_ramp, c(irow,:));
    % velocity = [velocity; vel];
    % FR = [FR; meanMaxFR'];
    % stepGCFRplots(P_step);
    stairGCFRplots(P_stair);

    % chirpGCFRplots(P_chirp(1));
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