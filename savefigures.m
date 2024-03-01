% function savefigures(filename, stim_name, figurehandle, date)
function savefigures(P, addlabel, figureHandle, imageFlag, dir_path)
    
%     dir_path = 'F:\Chinmayee\Work\Recordings\Intracellular recording\analysis';
   
    cd(dir_path)
%     svg_name = sprintf("%s_%s_%s", string(join(split(date,'.'),'_')), filename, stim_name);
%     prompt = 'Enter file name';
%     svg_name = input(prompt);
    % title(gca, "")
    if imageFlag == "svg" 
        ext = '.pdf';
        imageType = "vector";
    
    else
        ext = '.png';
        imageType = "image";
    end
    
    filename = join([string(P.date) "_" P.filename "_" P.stim_name "_" addlabel ext], "");
    % filename = join([string(P.date) "_" P.filename "_" addlabel ext], "");
%     imageType = 'image'
    if imageFlag == "fig"
        ext = ".fig";
        filename = join([string(P.date) "_" P.filename "_" P.stim_name "_" addlabel ext], "");
        % filename = join([string(P.date) "_" P.filename "_" addlabel ext], "");
        savefig(figureHandle, filename, "compact");
        return
    end
    exportgraphics(figureHandle, filename, 'ContentType', imageType, 'Resolution',150);
%     fig_name = sprintf("%s_%s.fig", filename, stim_name);
    
%     savefig(figurehandle, fig_name);
%     saveas(figurehandle, svg_name, 'svg');

end