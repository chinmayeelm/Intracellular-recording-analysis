% function savefigures(filename, stim_name, figurehandle, date)
function savefigures(figurehandle)
    
%     dir_path = 'F:\Chinmayee\Work\Recordings\Intracellular recording\analysis';
    dir_path = 'E:\Recordings\analysis\figure for report';
    cd(dir_path)
%     svg_name = sprintf("%s_%s_%s", string(join(split(date,'.'),'_')), filename, stim_name);
    prompt = 'Enter file name';
    svg_name = input(prompt);
    
%     fig_name = sprintf("%s_%s.fig", filename, stim_name);
    
%     savefig(figurehandle, fig_name);
    saveas(figurehandle, svg_name, 'svg');

end