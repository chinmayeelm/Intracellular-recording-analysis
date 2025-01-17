%% Quantification of adaptation
% 250 ms ramp protocols
% GCFR from peak + 3s 
% power2 curve fit
% rsquare >= 0.8 accepted

clear;
close all;
cd 'D:\Work\Code\Intracellular-recording-analysis\LUTs';
% rampData = readlines('velocityEncodingNeurons_2.txt');
rampData = readlines('positionEncodingNeurons.txt');
rampData = replace(rampData, "step", "ramp");
% rampData = velEncodingNeurons;
% rampData = posEncodingNeurons;
dataDirectory_ramp = extractBefore(rampData, '_');
filename_ramp = extractAfter(rampData, "_");
nfiles = length(filename_ramp);
rsquares = [];
adaptation_coeffs=[];
gcfrSSdur = 3;

for irow = 1:nfiles
    expt_date = replace(dataDirectory_ramp(irow), '-','.');
    filenameParts = split(filename_ramp(irow), '_');
    mothId = filenameParts(1);
    if exist(join(['D:\Work\Recordings\' string(expt_date) '\' 'raw\' mothId '\'], ''), 'file') ~=0
        P_ramp = getStructP(dataDirectory_ramp(irow), filename_ramp(irow),0);

        for i = 1:length(P_ramp)

            stim_name = split(P_ramp(i).stim_name);
            delta_t = str2double(stim_name(2));

            if delta_t == 0.2500
                
                t = P_ramp(i).time(1:P_ramp(i).single_trial_length);
                [~,maxFRLoc] = max(P_ramp(i).avg_gcfr);
                stim = P_ramp(i).mean_movement(maxFRLoc:maxFRLoc+3*P_ramp(i).fs);
                gcfr = P_ramp(i).avg_gcfr(maxFRLoc:maxFRLoc+3*P_ramp(i).fs); 
                t1 = linspace((1/P_ramp(i).fs), 3, length(gcfr));
                

                [xData, yData] = prepareCurveData( t1', gcfr );
                [fitresult, gof] = fit( xData, yData, 'power1' );
                adaptation_coeff = fitresult.b
                rsquare = gof.rsquare

                if rsquare >=0.8
                    rsquares = [rsquares rsquare];
                    adaptation_coeffs=[adaptation_coeffs adaptation_coeff];
                    % figure;
                    % subplot(2,1,1);
                    % plot(t, P_ramp(i).mean_movement,"Color",[0 0 0 0.5], LineWidth=1); hold on;
                    % xline(t(onLoc), 'g--');
                    % xline(t(offLoc), 'r--');
                    % ylabel("Position (deg)");
                    % xlabel('Time (s)');
                    % 
                    % subplot(2,1,2);
                    % plot(t,P_ramp(i).avg_gcfr,"Color",[0.8500 0.3250 0.0980 0.3], "LineWidth",1); hold on;
                    % xline(t(onLoc), 'g--');
                    % xline(t(offLoc), 'r--');
                    % ylabel('Normalised firing rate');
                    % xlabel('Time (s)');
                end
                

                % figure();
                % plot( fitresult, xData, yData ); hold on;
                % str = sprintf(" k = %0.3f \n rsquare = %0.3f", adaptation_coeff, rsquare);
                % text(t1(end)-2,max(gcfr)-10, str);
                % title(replace([P_ramp(i).date P_ramp(i).filename],'_',' '));
                % % ylim([0 1]);
                % % xlim([0 3.5]);
                % ylabel('Mean firing rate');
                % xlabel('Time (s)');
                % legend off;

            end
        end
    end

end

%%


figure; histogram(rsquares,'Normalization','probability','NumBins',20,'FaceColor','#0062BD', 'FaceAlpha', 0.2, 'EdgeAlpha',0.5)
xlim([0 1]);
ylabel('Probability');
xlabel('rsquares');


figure; histogram(adaptation_coeffs,'Normalization','probability','NumBins',20,'FaceColor','#0062BD', 'FaceAlpha', 0.2, 'EdgeAlpha',0.5)
ylabel('Probability');
xlabel('Adaptation Coefficient');
box off

figure; boxplot(adaptation_coeffs);
hold on;
swarmchart(ones([1,length(adaptation_coeffs)]), adaptation_coeffs, 'o', "filled", "XJitterWidth", 0.1,"MarkerFaceColor",'#D95319')
xticklabels('All neurons');
% ylim([-1.5 0]);
% yticks(-1.5:0.25:0);
ylabel('Adaptation coefficient');
box off;