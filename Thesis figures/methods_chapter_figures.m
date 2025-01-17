    %% Methods Chapter

axisLineWidth = 1;
axisFont = 'Calibri';
labelFont = 'Calibri';
tickLabelSize = 10;
labelFontSize = 12;


%% Filtering neural activity and spike detection
% 21.06.2022 M2_N4_ramp

figure();
% raw data
ax1 = subplot(2,1,1); plot(time(1:fs), 100*rec_data(1:fs), 'k', 'LineWidth',1);
ax1.FontName = 'Calibri';
ax1.FontSize = 8;
ax1.LineWidth = 1;
ylabel({ 'Raw', 'Membrane', 'potential', '(mV)'}, 'Rotation',0, 'FontSize', 12);
ylim([-80 0]);
% xlim('padded');
[hy,hx] = offsetaxis(ax1, 'y', 0.05, 'x', 0.1, 'yloc', 'l', 'xloc', 'b');
hx.Visible = 'off';
box off;
% ax1_pos = get(ax1, 'Position');
% ax1.Visible = 'off';
% ax = axes;
% ax.XAxis.Position = ax1_pos(1) - 0.02;
% ax.YAxis.Position = ax1_pos(2) - 0.02;
% set(ax1, 'Position', [ax1_pos(1)+0.02 ax1_pos(2) ax1_pos(3) ax1_pos(4)]);

% filtered data and Spike detection
ax2 = subplot(2,1,2); plot(time(1:fs), 100*filtered_data_bp(1:fs), 'k', 'LineWidth',1); hold on;
[p,l] = findpeaks(100*filtered_data_bp(1:fs), "MinPeakHeight",0.3*max(100*filtered_data_bp(1:fs)), "MinPeakDistance", 0.003*fs);
plot(l/fs, p, '.', 'MarkerEdgeColor', 'r');
ax2.FontName = 'Calibri';
ax2.FontSize = 8;
ax2.LineWidth = 1;
yline(0.25*max(100*filtered_data_bp(1:fs)), 'b--', 'LineWidth',1);
ylabel({'Filtered', 'Membrane', 'potential', '(mV)'}, 'Rotation',0, 'FontSize', 12);
xlabel('Time (s)', 'FontSize',12);
box off;
% xlim('padded');
[hy,hx] = offsetaxis(ax2, 'y', 0.05, 'x', 0.1, 'yloc', 'l', 'xloc', 'b');
% ax2_pos_x = get(hx, 'Position');
% set(hx, 'Position', [ax2_pos_x(1)+0.02 ax2_pos_x(2)-0.02 ax2_pos_x(3) ax2_pos_x(4)]);
% ax2_pos_y = get(hy, 'Position');
% set(hy, 'Position', [ax2_pos_y(1)+0.02 ax2_pos_y(2) ax2_pos_y(3) ax2_pos_y(4)]);

% linkaxes([ax1,ax2], 'x');


%% Generated stimuli

%12.07.2022 M1_N3


dataDirectory = '2022.07.20';
file_names = ["stair"; "ramp"; "step"; "chirp"; "blwgn"; "sweep_sqr"; "sine"];
protocol_names = ["stair"; "ramp"; "step"; "chirp"; "blwgn"; "sweep"; "sin"; "sqr"];
nFiles = length(file_names);
nPrtcls = length(protocol_names);
figure;

for iFile = 1:nFiles
    filename = sprintf('M1_N3_%s', file_names(iFile));
    disp(filename);
    P = getStructP(dataDirectory, filename,0,1);

    for iPrtcl = 1:nPrtcls
        for i=1:length(P)
            if contains(P(i).stim_name, protocol_names(iPrtcl))
              
                subplot(nPrtcls,1,iPrtcl);
                plot(P(i).time(1:P(i).single_trial_length), P(i).intendedStimulus(1,:), 'k'); hold on;
                line([1 2], [0 0],'Color', 'k', 'LineWidth', 2);
                text(4.5, 0.002, '1 s');
                ax = gca;
                ax.XAxis.Visible = "off";
                ax.YAxis.Visible = "off";
                ax.YAxis.Limits = [-0.025 0.025];
                box off;
            end
        end
    end
end

%% Comparison of HES and generated signal

%12.07.2022 M1_N3


dataDirectory = '2022.07.20';
file_names = ["stair"; "ramp"; "step"; "sweep_sqr"; "sine"];
% file_names = "sweep_sqr";
protocol_names = ["stair"; "ramp"; "step"; "sweep"; "sqr"; "sin"];

nFiles = length(file_names);
nPrtcls = length(protocol_names);


for iFile = 1:nFiles
    filename = sprintf('M1_N3_%s', file_names(iFile));
    disp(filename);
    P = getStructP(dataDirectory, filename,[nan nan],1);
    figure;
    for iPrtcl = 1:nPrtcls
        for i=1:length(P)
            if contains(P(i).stim_name, protocol_names(iPrtcl))
              
                subplot(nPrtcls,1,iPrtcl);
                scalingFactor = max(P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs))) - min(P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs)));
                intended_pp = max(P(i).intendedStimulus(1,1*P(i).fs:end-(1*P(i).fs))) - min(P(i).intendedStimulus(1,1*P(i).fs:end-(1*P(i).fs)));
                % yyaxis right; 
                plot(P(i).time(1*P(i).fs:P(i).single_trial_length-(1*P(i).fs)), -(scalingFactor/intended_pp)*P(i).intendedStimulus(1,1*P(i).fs:end-(1*P(i).fs)), 'k', 'LineStyle','-','Marker','none'); hold on; 
                set(gca, "YLimMode", "auto");    
                % yyaxis left; 
                plot(P(i).time(1*P(i).fs:P(i).single_trial_length-(1*P(i).fs)), P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs)), 'Color', [0.6, 0.2,0], 'LineStyle','-','Marker','none'); 
                % sdfill(P(i).time(1*P(i).fs:P(i).single_trial_length-(1*P(i).fs)), P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs)), std(P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs)),[],1), [0.6, 0.2,0])
                % set(gca, "YLimMode", "auto");
                
                % ylim([min(P(i).antennal_movement(1,1*P(i).fs:end - (1*P(i).fs))) max(P(i).antennal_movement(1,1*P(i).fs:end - (1*P(i).fs)))]);
                line([1 5], [0 0],'Color', 'k', 'LineWidth', 2);
                line([1 1],[0 0.5], 'Color', 'k', 'LineWidth', 2); 
                text(2.5, -0.001, '4 s');
                text(1,0.002, '0.5deg', 'Interpreter','tex')
                ax = gca;
                ax.XAxis.Visible = "off";
                ax.YAxis(1).Visible = "off";
                % ax.YAxis(2).Visible = "off";
                % ax.YAxis.Limits = [-0.025 0.025];
                box off;
            end
        end
    end
end

%% Chirp/blwgn zoom comparison with generate

dataDirectory = "2022.07.20";
filename = "M1_N3_chirp";

P = getStructP(dataDirectory, filename,[nan nan],1);
% plot_data(P(1));

i=1%:length(P)
figure;

scalingFactor = max(P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs))) - min(P(i).antennal_movement(1,1*P(i).fs:end-(1*P(i).fs)));
intended_pp = max(P(i).intendedStimulus(1,1*P(i).fs:end-(1*P(i).fs))) - min(P(i).intendedStimulus(1,1*P(i).fs:end-(1*P(i).fs)));

% subplot(2,1,1);
subplot(3,1,1);
hold on;
patch([P(i).OFF_dur P(i).OFF_dur P(i).OFF_dur+P(i).ON_dur P(i).OFF_dur+P(i).ON_dur], [-2 2 2 -2],[0 0 0], 'FaceAlpha' , 0.1, 'EdgeColor', 'none');
plot(P(i).time(1:P(i).single_trial_length), (scalingFactor/intended_pp)*P(i).intendedStimulus, 'Color', 'k');
for itrial = 1%:P(i).complete_trials
    plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement(itrial, :), 'Color', [0.6, 0.2,0], 'LineWidth',1);
end


% line([4 5],[-0.6 -0.6], 'Color','k', 'LineWidth', 2);
% text(4, -0.75, '1 s');
% line([4 4],[-0.6 -0.1], 'Color','k', 'LineWidth', 2);
% text(3, -0.3, '0.5deg' );


xlim_main = [0 2*P(i).OFF_dur+P(i).ON_dur];
ylim([-2 2]);
xlim(xlim_main);
ax = gca;
ax.FontName = 'Calibri'; 
% ylabel('Angular position (deg)');
box off;

% subplot(2,1,2);
subplot(3,1,2);

hold on;
patch([P(i).OFF_dur P(i).OFF_dur P(i).OFF_dur+P(i).ON_dur P(i).OFF_dur+P(i).ON_dur], [-2 2 2 -2],[0 0 0], 'FaceAlpha' , 0.1, 'EdgeColor', 'none');
plot(P(i).time(1:P(i).single_trial_length), (scalingFactor/intended_pp)*P(i).intendedStimulus, 'Color','k');
for itrial = 1%:5%P(i).complete_trials
    plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement(itrial, :), 'Color', [0.6, 0.2,0], 'LineWidth',1);
end

% line([4.5 5],[-0.6 -0.6], 'Color','k', 'LineWidth', 2);
% text(4.5, -0.75, '0.5 s');
% line([4.5 4.5],[-0.6 -0.1], 'Color','k', 'LineWidth', 2);
% text(4.25, -0.4, '0.5deg' );
ylim([-2 2]);
xlim_zoom = miscFuncs.range_converter(xlim_main, P(i).OFF_dur, 100);
xlim(xlim_zoom);
xlabel('Time (s)');
ylabel('Angular position (deg)');
ax = gca;
ax.FontName = 'Calibri'; 
box off;

subplot(3,1,3);
hold on;
patch([5 5 20 20], [-1 1 1 -1],[0 0 0], 'FaceAlpha' , 0.1, 'EdgeColor', 'none');
plot(P(i).time(1:P(i).single_trial_length), (scalingFactor/intended_pp)*P(i).intendedStimulus, 'Color', 'k');
for itrial = 1%:P(i).complete_trials
    plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement(itrial, :), 'Color', [0.6, 0.2,0], 'LineWidth',1);
end

% line([4.5 5],[-0.6 -0.6], 'Color','k', 'LineWidth', 2);
% text(4.5, -0.75, '0.5 s');
% line([4.5 4.5],[-0.6 -0.1], 'Color','k', 'LineWidth', 2);
% text(4.25, -0.4, '0.5deg' );
ylim([-1 1]);
xlim_zoom = miscFuncs.range_converter(xlim_main, P(i).OFF_dur+P(i).ON_dur, 10);
xlim(xlim_zoom);
xlabel('Time (s)');
% ylabel('Angular position (deg)');
ax = gca;
ax.FontName = 'Calibri'; 
box off;


%% Comparison of frequency spectrum WN
dataDirectory = "2022.07.20";
filename = "M1_N3_blwgn";

P = getStructP(dataDirectory, filename,[nan nan],1);

startPt = P.OFF_dur*P.fs+1;
stopPt = (P.OFF_dur+P.ON_dur)*P.fs;
stim_length = length(P.intendedStimulus(startPt:stopPt));
win_length = stim_length/10;
overlap = 0.9*win_length;
[pxx_gen, f_gen] = pwelch(P.intendedStimulus(1,startPt:stopPt),win_length,overlap,[], P.fs);
[pxx_hes, f_hes] = pwelch(P.antennal_movement(1,startPt:stopPt),win_length,overlap,[], P.fs);

% [pxx_gen, f_gen] = pwelch(P.intendedStimulus(1,P.OFF_dur*P.fs+1:(P.OFF_dur+P.ON_dur)*P.fs),[],[],[], P.fs);
% [pxx_hes, f_hes] = pwelch(P.antennal_movement(1,P.OFF_dur*P.fs+1:(P.OFF_dur+P.ON_dur)*P.fs),[],[],[], P.fs);

figure;
hold on;
plot(f_gen, 10*log10(pxx_gen/max(pxx_gen)), 'k');
plot(f_hes, 10*log10(pxx_hes/max(pxx_hes)), 'Color', [0.6, 0.2,0]);
ylabel('Normalized Power spectral density (dB/Hz)');
xlabel('Frequency (Hz)');
xlim([0 400]);
set(gca, 'FontName', 'Calibri');

%% Comparison of frequency spectrum WN
dataDirectory = "2022.07.20";
filename = "M1_N3_chirp";

P = getStructP(dataDirectory, filename,[nan nan],1);
P = P(2);
startPt = P.OFF_dur*P.fs+1;
stopPt = (P.OFF_dur+P.ON_dur)*P.fs;
stim_length = length(P.intendedStimulus(startPt:stopPt));
win_length = stim_length/10;
overlap = 0.8*win_length;
[pxx_gen, f_gen] = pwelch(P.intendedStimulus(1,startPt:stopPt),win_length,overlap,[], P.fs);
[pxx_hes, f_hes] = pwelch(P.antennal_movement(1,startPt:stopPt),win_length,overlap,[], P.fs);

figure;
hold on;
plot(f_gen, 10*log10(pxx_gen/max(pxx_gen)), 'k');
plot(f_hes, 10*log10(pxx_hes/max(pxx_hes)), 'Color', [0.6, 0.2,0]);
ylabel('Normalized Power spectral density (dB/Hz)');
xlabel('Frequency (Hz)');
xlim([0 200]);
set(gca, 'FontName', 'Calibri');

%% Repeatability

dataDirectory = '2022.07.20';
filename = "M1_N3_blwgn"
P = getStructP(dataDirectory, filename,[nan nan],1);

    
i=1%:length(P)
figure;

subplot(2,1,1);
hold on;
patch([P(i).OFF_dur P(i).OFF_dur P(i).OFF_dur+P(i).ON_dur P(i).OFF_dur+P(i).ON_dur], [-2 2 2 -2],[0 0 0], 'FaceAlpha' , 0.1, 'EdgeColor', 'none');
for itrial = 1:5%P(i).complete_trials
    plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement(itrial, :), 'Color', [0.6, 0.2,0, 0.2], 'LineWidth',1);
end
plot(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, 'Color', [0.6, 0.2,0], 'LineWidth',1);

% line([4 5],[-0.6 -0.6], 'Color','k', 'LineWidth', 2);
% text(4, -0.75, '1 s');
% line([4 4],[-0.6 -0.1], 'Color','k', 'LineWidth', 2);
% text(3, -0.3, '0.5deg' );


xlim_main = [0 2*P(i).OFF_dur+P(i).ON_dur];

xlim(xlim_main);
ax = gca;
ax.FontName = 'Calibri'; 
ylabel('Angular position (deg)');
box off;


subplot(2,1,2);

hold on;
patch([P(i).OFF_dur P(i).OFF_dur P(i).OFF_dur+P(i).ON_dur P(i).OFF_dur+P(i).ON_dur], [-2 2 2 -2],[0 0 0], 'FaceAlpha' , 0.1, 'EdgeColor', 'none');
for itrial = 1:5%P(i).complete_trials
    plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement(itrial, :), 'Color', [0.6, 0.2,0, 0.2], 'LineWidth',1);
end
plot(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, 'Color', [0.6, 0.2,0], 'LineWidth',1);
% line([4.5 5],[-0.6 -0.6], 'Color','k', 'LineWidth', 2);
% text(4.5, -0.75, '0.5 s');
% line([4.5 4.5],[-0.6 -0.1], 'Color','k', 'LineWidth', 2);
% text(4.25, -0.4, '0.5deg' );
ylim([-1 1]);
xlim_zoom_start = 4.9;
% xlim_zoom = [xlim_zoom_start xlim_zoom_start + (P(i).OFF_dur-xlim_zoom_start)*diff(xlim_main)/(P(i).OFF_dur - xlim_main(1))]
xlim_zoom = miscFuncs.range_converter(xlim_main, 5, 30);
xlim(xlim_zoom);
xlabel('Time (s)');
ylabel('Angular position (deg)');
ax = gca;
ax.FontName = 'Calibri'; 
box off;

% subplot(3,1,3);
% 
% hold on;
% patch([5 5 20 20], [-1 1 1 -1],[0 0 0], 'FaceAlpha' , 0.2, 'EdgeColor', 'none');
% for itrial = 1:P(i).complete_trials
%     plot(P(i).time(1:P(i).single_trial_length), P(i).antennal_movement(itrial, :), 'Color', [0.6, 0.2,0, 0.3], 'LineWidth',1);
% end
% plot(P(i).time(1:P(i).single_trial_length), P(i).mean_movement, 'Color', [0.6, 0.2,0], 'LineWidth',1);
% % line([4.5 5],[-0.6 -0.6], 'Color','k', 'LineWidth', 2);
% % text(4.5, -0.75, '0.5 s');
% % line([4.5 4.5],[-0.6 -0.1], 'Color','k', 'LineWidth', 2);
% % text(4.25, -0.4, '0.5deg' );
% ylim([-1 1]);
% xlim_zoom_start = 19.5;
% xlim_zoom = [19.6 20.1];
% xlim(xlim_zoom);
% xlabel('Time (s)');
% ylabel('Angular position (deg)');
% ax = gca;
% ax.FontName = 'Calibri'; 
% box off;
        


%% Sensitivity and SNR
% Need another reference for such analysis




%% STD across trials
% 22.07.2022

dataDirectory = "2022.06.21";
filename = "M2_N4_blwgn";
P = getStructP(dataDirectory, filename,[nan nan],1);

figure;
plot(P(1).time(1:P(1).single_trial_length), P(1).intendedStimulus(1,:), 'LineWidth',1); hold on;
yyaxis right; plot(P(1).time(1:P(1).single_trial_length), P(1).antennal_movement(1,:)); 

figure;
plot(P(1).time(1:P(1).single_trial_length), [P(1).antennal_movement]); 

figure;
plot(P(1).time(1:P(1).single_trial_length), [P(1).antennal_movement(1:2,:)]); hold on;
% yyaxis right; plot(P(1).time(1:P(1).single_trial_length), P(1).intendedStimulus(1,:));

figure;
% fft_stim_(P(1).intendedStimulus(1,:), P(1).fs);
fft_stim_(P.antennal_movement(1,:), P.fs); hold on;
fft_stim_(P.antennal_movement(2,:), P.fs); 


%%
figure;
for irow = 1:height(LUT_blwgn_fs_20k)
  
    dataDirectory = LUT_blwgn_fs_20k.expt_date(irow);
    filename = LUT_blwgn_fs_20k.filename(irow);
    
    P = getStructP(dataDirectory, filename, [nan nan],1);
    t = P(1).time(1:P(1).single_trial_length);
    plot(t, P(1).intendedStimulus(1,:), 'DisplayName',replace(join([dataDirectory filename], " "),"_", " "));
    hold on;

end


%%
dataDirectory = "2022.06.21";
filename = "M2_N4_blwgn";
P = getStructP(dataDirectory, filename,[nan nan],1);
startPt = P(1).OFF_dur*P(1).fs+1;
stopPt = (P(1).OFF_dur + P(1).ON_dur)*P(1).fs;

 for itrial = 1:P.complete_trials
        
    stim_hes_unfilt = P.hes_data_unfilt(itrial,:) - P.hes_data_unfilt(itrial,1);
    [x,y] = butter(10, 300/(P(1).fs/2), "low");
    hes_filt = filtfilt(x,y, stim_hes_unfilt);
    stim_hes_filt = hes_filt(startPt : stopPt);

    stim_ifb = P.stim_ifb(itrial,startPt:stopPt) - P.stim_ifb(itrial,1);
    % stim_intended = P.intendedStimulus(itrial,startPt:stopPt);
    gcfr = P.gcfr(itrial, startPt:stopPt);
    c_ifb_gcfr = corrcoef(stim_ifb, gcfr);
    c_ifb = c_ifb_gcfr(2)
    c_hes_gcfr = corrcoef(stim_hes_filt, gcfr);
    c_hes = c_hes_gcfr(2)

end



