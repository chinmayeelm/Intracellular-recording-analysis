% function savefigures(filename, stim_name, figurehandle, date)
function savefigures(P, addlabel, figureHandle, imageFlag)
    
%     dir_path = 'F:\Chinmayee\Work\Recordings\Intracellular recording\analysis';
    dir_path = 'F:\Work\Figures for presentation';
    cd(dir_path)
%     svg_name = sprintf("%s_%s_%s", string(join(split(date,'.'),'_')), filename, stim_name);
%     prompt = 'Enter file name';
%     svg_name = input(prompt);
    title(gca, "")
    if imageFlag == 0
        ext = '.pdf';
        imageType = "vector";
    else
        ext = '.png';
        imageType = "image";
    end
    
    filename = join([string(P(1).date) "_" P(1).filename "_" addlabel ext], "")
%     imageType = 'image'
    exportgraphics(figureHandle, filename, 'ContentType', imageType, 'Resolution',150);
%     fig_name = sprintf("%s_%s.fig", filename, stim_name);
    
%     savefig(figurehandle, fig_name);
%     saveas(figurehandle, svg_name, 'svg');

end