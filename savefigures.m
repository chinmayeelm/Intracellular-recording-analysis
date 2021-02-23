function savefigures(filename, stim_name, figurehandle, date)
    
%     dir_path = 'F:\Chinmayee\Work\Recordings\Intracellular recording\analysis';
    
    svg_name = sprintf("%s_%s_%s", string(join(split(date,'.'),'_')), filename, stim_name);
    
%     fig_name = sprintf("%s_%s.fig", filename, stim_name);
    
%     savefig(figurehandle, fig_name);
    saveas(figurehandle, svg_name, 'fig');

end