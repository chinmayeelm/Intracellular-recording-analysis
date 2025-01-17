LUT_path = 'D:\Work\Code\Intracellular-recording-analysis\LUTs\JOrecLUT-01062022onwards.mat';
load(LUT_path);


for irow=81%:height(LUT_intra)
    irow
    close all;
    dataDirectory = LUT_intra.expt_date(irow)
    filename = LUT_intra.filename(irow)




    % if contains(filename, "wgn")
    %     continue;
    % end

    % if contains(filename, "ramp")
            
    LUTrow = LUT_intra((LUT_intra.expt_date==dataDirectory & LUT_intra.filename == filename),:);
    P = getStructP_LUT(dataDirectory, filename, LUTrow, [LUT_intra.start_clip(irow) LUT_intra.stop_clip(irow)],1);
    protocolPlot(P);


    % else
    %     continue;
    % end
    % prompt = "Should I plot protocols?";
    % plot_flag = input(prompt);
    % 
    % if plot_flag == 1
    %     protocolPlot(P);
    % else
    %     continue;
    % end

    % if contains(filename, "ramp")
    %     protocolPlot(P);
    %     % rampGCFRplots(P);
    % % elseif contains(filename, "step")
    % %     stepGCFRplots(P);
    % end

    pause;


end


