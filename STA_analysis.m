function STA = STA_analysis(stimulus, raster, stim_window, fs, varargin)
    
    narginchk(4, 6);

    if nargin==6
        STE = varargin{end-1};
        PSE = varargin{end};
    else
        % STE = getIsolatedSTE(stimulus, raster, stim_window, fs);
        % PSE = getPSE(stimulus, stim_window, fs, size(STE,1)*100);
    end

    

    avg_stim = mean(PSE,1);
    STA = mean(STE,1)-avg_stim;
    % STA = mean(STE,1);
    
    %{
    % sd = std(STE-avg_stim,[],1);
    sd = std(STE,[],1);
    % sem = sd./sqrt(size(STE,1));
    t_STA = linspace(-(stim_window*1000),0,length(STA));%-100:0.1:0;

    figure(); 
    plot(t_STA, STA); hold on;
    sdfill(t_STA, STA, sd, 'k');

    % subtitle ('Spike triggered average');
    ylabel 'Antennal movement (deg)';
    xlabel 'Time (ms)';    
    box off;
    %}

end




