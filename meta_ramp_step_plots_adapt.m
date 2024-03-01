cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
%
% filepath = readlines('stepDataList.txt');
% filepath = readlines('allStairList.txt');
% filepath = readlines('chirpPartialList.txt');
% filepath = readlines('blwgnValidData.txt');
% filepath = readlines('ramp_forAdaptation.txt');
% filepath = readlines('validRampList.txt');
% filepath = readlines('rampDataListNew.txt');

dataDirectory = extractBefore(filepath, '_');
filename = extractAfter(filepath, "_");
nfiles = length(filename);

% c = parula(nfiles);
% T_tonic = table(); %('neuronID', 'velocity', 'Tphasic', 'Ttonic', 'pvalTphasic', 'pvalTtonic');
% T_ramp_gcfr_std = table();
T= table();


for irow = 1:nfiles
    irow
    close all;

    P = getStructP(dataDirectory(irow), filename(irow), [nan nan], 1);


    % 
    

    % pause;



    % figure('WindowState','fullscreen'); protocolPlot(P);
    % savefigures(P(1), "traces", gcf, "fig", 'D:\Work\Figures for presentation\step_traces');
    % savefigures(P(1), "traces", gcf, "png", 'D:\Work\Figures for presentation\step_traces');

    %Estimating adaptation

    % if contains(P(1).stim_name, "ramp")
    % stim_name = string(extractfield(P, 'stim_name'));
    % ramp_dur = str2double(extractAfter(stim_name, "ramp "));
    % [~,idx] = sort(ramp_dur);
    % P = P(idx);
    % protocolPlot(P);
    % baseline_sd = [];
    % ss_sd = [];
    % for i=1:length(P)
    %     [onLoc, offLoc] = miscFuncs.findSSbounds(P(i).mean_movement, 0.9, 10, P(i).fs);
    %     ssBounds = [offLoc-1.5*P(i).fs  offLoc-0.5*P(i).fs];
    %     baselineBounds = [(P(i).OFF_dur-1.5)*P(i).fs+1 (P(i).OFF_dur-0.5)*P(i).fs];
    % 
    %     baseline_sd = [baseline_sd; std(P(i).gcfr(:, baselineBounds(1):baselineBounds(2)), [], "all")];
    %     ss_sd = [ss_sd; std(P(i).gcfr(:, ssBounds(1):ssBounds(2)), [], "all")];
    % 
    % 
    % end
    %     neuronName = extractBefore(P(i).filename, "_ramp");
    %     nID = replace(join([P(i).date neuronName], ""), "_", "");
    %     nPrtcls = length(baseline_sd);
    %     neuronID = repmat(nID, nPrtcls,1);
    % 
    %     Trow = table(baseline_sd, ss_sd, neuronID);
    %     T_ramp_gcfr_std = [T_ramp_gcfr_std; Trow];
    %{
        Tphasic = [];
        rsqExp = [];
        velocity = [];
        peakFR = [];

        for i=1:length(P)

            [max_FR, vel, fitresult, gof ] = estimatePhasicAdaptation(P(i));
            % Tphasic = [Tphasic; (1/mdl_phasic_exp.Coefficients.Estimate(2))];
            % pvalTphasic = [pvalTphasic; mdl_phasic_exp.Coefficients.pValue(2)];
            % rsqExp = [rsqExp; mdl_phasic_exp.Rsquared.Ordinary];
            Tphasic = [Tphasic; (1/fitresult.b)];
            rsqExp = [rsqExp; gof.rsquare];
            
            velocity = [velocity; vel];
            peakFR = [peakFR; max_FR];

        end

        neuronName = extractBefore(P(i).filename, "_ramp");
        nID = replace(join([P(i).date neuronName], ""), "_", "");
        nPrtcls = length(velocity);
        neuronID = repmat(nID, nPrtcls,1);
        Trow = table(neuronID, velocity, Tphasic, peakFR, rsqExp);
        T_phasic = [T_phasic;Trow];
    end

    if contains(P(1).stim_name, "step")

        stim_name = string(extractfield(P, 'stim_name'));
        px = str2double(extractAfter(stim_name, "amp_"));
        [~,idx] = sort(px);
        P = P(idx);

        Ttonic = [];
        pvalTtonic = [];
        rsqExp = [];
        position = [];
        ssFR = [];
        

        for i=1:length(P)
            [ss_FR, pos, fitresult, gof ] = estimateTonicAdaptation(P(i));
            % Ttonic = [Ttonic; (1/mdl_tonic_exp.Coefficients.Estimate(2))];
            % pvalTtonic = [pvalTtonic; mdl_tonic_exp.Coefficients.pValue(2)];
            % rsqExp = [rsqExp; mdl_tonic_exp.Rsquared.Ordinary];
            Ttonic = [Ttonic; (1/fitresult.b)];
            rsqExp = [rsqExp; gof.rsquare];

            position = [position; pos];
            ssFR = [ssFR; ss_FR];

        end

        neuronName = extractBefore(P(i).filename, "_step");

        nID = replace(join([P(i).date neuronName], ""), "_", "");
        nPrtcls = length(position);
        neuronID = repmat(nID, nPrtcls,1);
        Trow = table(neuronID, position, Ttonic, ssFR, rsqExp);
        T_tonic = [T_tonic;Trow];
    end
    %}

    % [gcfrDecay, tDecay, phasicGCFR, tPhasic] = getPhasicDecayGCFR(P(1));
    % [dist1,dist2] = distanceToFit(t1, gcfr, rampEndIdx/P(i).fs, P(i).fs);
    % [velocity, amp, maxFR] = rampGCFRplots(P,c(irow, :));

    % bounds = [(P(1).OFF_dur-1)*P(1).fs (P(1).OFF_dur+5)*P(1).fs];
    % if contains(P(1).stim_name, "1")
    % % Trow = table({velocity}, {maxFR});
    % Trow = table(P(1).mean_movement(bounds(1):bounds(2)),P(1).avg_gcfr(bounds(1):bounds(2)));
    % T(irow, :) = Trow;
    % end


    % [pvalA, pvalAB, steadystateFR, baselineFR_mat] = stepGCFRplots(P,c(irow, :));

    % stairGCFRplots(P(i));
    % if P(end).max_chirp_frq == 150 && P(end).ON_dur == 15 && P(end).OFF_dur == 5
    %     irow
    %     chirpGCFRplots(P(end));
    % end
    % figure;



    % pause;

end

