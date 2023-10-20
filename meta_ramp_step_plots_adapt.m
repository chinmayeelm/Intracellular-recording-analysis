cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
% rampData = readlines('rampDataList.txt');
% dataDirectory_ramp = extractBefore(rampData, '_');
% filename_ramp = extractAfter(rampData, "_");
% nfiles = length(filename_ramp);

% stepData = readlines('stepDataList.txt');
% % stepData = stepData_(posEncodingNeurons,:);
% dataDirectory_step = extractBefore(stepData, '_');
% filename_step = extractAfter(stepData, "_");
% nfiles = length(filename_step);

stairData = readlines('allStairList.txt');
dataDirectory_stair = extractBefore(stairData, '_');
filename_stair = extractAfter(stairData, "_");
nfiles = length(filename_stair);

% chirpData = readlines('chirpPartialList.txt');
% chirpData = chirp_pos;
% dataDirectory_chirp = extractBefore(chirpData, '_');
% filename_chirp = extractAfter(chirpData, "_");
% nfiles = length(filename_chirp);

% wnData = readlines('blwgnValidData.txt');
% wnData = blwgn_vel;
% dataDirectory_wn = extractBefore(wnData, '_');
% filename_wn = extractAfter(wnData, "_");
% nfiles = length(filename_wn);

% adaptation_coeff = zeros(1, nfiles);
% rsquare = zeros(1, nfiles);

% ignorefiles = [];
% velocity = [];
% FR = [];
% velocityEncodingNeurons = [];
% staList = [];
% posEncodingNeurons = [];

% figure;
c = parula(nfiles);
ci = 0;
% goodfits = [];
% slopes = [];
% rsquares = [];

labelFontSize = 14;
tickLabelSize = 12;

for irow = 1:nfiles

    irow
    %     close all;

    ci = ci+1;
    % try
    % P_ramp = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),0);
    % P_step = getStructP(dataDirectory_step(irow), filename_step(irow),0);
    P_stair = getStructP(dataDirectory_stair(irow), filename_stair(irow),0);
    % P_chirp = getStructP(dataDirectory_chirp(irow), filename_chirp(irow),0);
    % me_chirp(irow),0);
    % P_wn = getStructP(dataDirectory_wn(irow), filename_wn(irow),0);
    % title(irow)

    %     S(irow) = jitter_func(P_wn, irow);
    % [vel, meanMaxFR, p1, rsq] = rampGCFRplots(P_ramp, c(irow,:));

    % if p1>= 30 && rsq >= 0.8
    %     velocityEncodingNeurons = [velocityEncodingNeurons rampData(irow)];
    % end
    % velocity = [velocity; vel];
    % FR = [FR; meanMaxFR'];
    % [pvalA, pvalAB, steadystateFR, baselineFR_mat] = stepGCFRplots(P_step);
    % 
    % if ~isempty(find(pvalA <= 0.008)) 
    %     pvalA;
    %     posEncodingNeurons = [posEncodingNeurons irow];
    % end
    % pause;
    stairGCFRplots(P_stair);
    cd('D:\Work\Recordings\plots\stair');
    saveas(gcf,num2str(irow),'png');
    close all;
    % chirpGCFRplots(P_chirp);
    % try
    % [adaptation_coeff(irow), rsquare(irow)] = calcAdaptation(P_ramp);
    % catch
    %     disp("Not enough duration to fit curve");
    %     ignorefiles = [ignorefiles irow];
    % end
    %
    % catch
    %     ignorefiles = [ignorefiles irow];
    % end
    % if P_wn.fs ~= 20000
    %     continue;
    % end
    % STA_window = 0.03;
    % start_stim = P_wn.OFF_dur*P_wn.fs;
    % stop_stim = (P_wn.ON_dur+P_wn.OFF_dur)*P_wn.fs;
    % STA  = STA_analysis(P_wn.raster, P_wn.antennal_movement, STA_window, P_wn.fs, start_stim, stop_stim);
    % % staList = [staList; STA];
    % t = linspace(-30,0, length(STA));
    % % figure;
    % plot(t, STA, 'LineWidth',1); hold on;
    % box off;
    % ylabel('Antennal movement (deg)');
    % xlabel('Time before spike (ms)');
    % title(replace(join([P_wn.date ' ' P_wn.filename], ''), '_', ' '));
    % pause;
    % close all;
    % end
    
    %{
    [xData, yData] = prepareCurveData( vel, meanMaxFR );
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );

    slopes = [slopes fitresult.p1];
    rsquares = [rsquares gof.rsquare];
    if fitresult.p1 >= 40 && gof.rsquare >= 0.9
        goodfits = [goodfits irow];
        alphaval = 1;
        c_alpha  = [0.4660 0.6740 0.1880 alphaval];
    else
        alphaval = 0.3;
        c_alpha  = [0 0 0 alphaval];
    end


    % plot(fitresult,'k--');  hold on;
    plot(vel, meanMaxFR, 'Color',c_alpha,'LineWidth',1); hold on;

%}
    % pause;
    
end
% ax = gca;
% ylabel('Mean peak firing rate (Hz)', 'FontSize',labelFontSize, 'FontName','Calibri');
% xlabel('Angular velocity (deg/s) (log10)', 'FontSize',labelFontSize, 'FontName','Calibri');
% ax.LineWidth =1;
% ax.FontSize = tickLabelSize;
% ax.FontName = 'Calibri';
% ax.Box = 'off';
