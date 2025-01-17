function [eVal, sig_evec, sig_evec_orth, sta,avg_stim, cov_matrix, stim_prior_cov, diff_cov] = cov_analysis(stimulus, raster, stim_window, fs, varargin)
    
    narginchk(4, 6);

    if nargin==6
        STE = varargin{end-1};
        PSE = varargin{end};
    else
        disp('Generating STE and PSE');
        STE = getIsolatedSTE(stimulus, raster, stim_window, fs);
        PSE = getPSE(stimulus, stim_window, fs, size(STE,1)*100);
    end
    
    avg_stim = mean(PSE, 1);
    stim_prior_cov = cov(PSE - avg_stim);
    STA = mean(STE - avg_stim, 1);
    sta = STA';
    nTrials = size(stimulus, 1);
    nspikes = size(STE, 1);
    cov_matrix = cov(STE - STA);


    % Null distribution of eigenvalues
    pattern_length = stim_window*fs;
    rand_stim = zeros(nspikes, pattern_length);
    nRepeats = 500;
    Ds_null = nan(nRepeats, pattern_length);

    for irepeat = 1:nRepeats
        rr = randi([1, size(stimulus, 2) - pattern_length], 1, nspikes);
        parfor j = 1:nspikes
            rand_stim(j, :) = stimulus(randperm(nTrials, 1), rr(j):rr(j) + pattern_length - 1);
        end
        random_cov = cov(rand_stim - STA);
        [~, D] = eig(random_cov - stim_prior_cov);
        Ds_null(irepeat,:) = diag(D)';
    end
    

    min_null = min(Ds_null, [], 'all')
    max_null = max(Ds_null, [], 'all')

    diff_cov = cov_matrix - stim_prior_cov;

    
    [V, D] = eig(diff_cov);
    [~, ind] = sort(abs(real(diag(D))), 'descend');
   
    Ds = D(ind, ind);
    Vs = V(:, ind);

    eVal = real(diag(Ds));
    eVal_sig_ind = find(eVal < min_null | eVal > max_null);
    % eVal_sig_ind = find(eVal > max_null);
    eVal_sig = eVal(eVal_sig_ind);
    sig_evec = eVal_sig' .* Vs(:, eVal_sig_ind);

    sig_evec_orth = sig_evec - ((dot(sig_evec, repmat(STA',1,size(sig_evec,2))) ./ norm(STA)^2) .* STA');


    % figure;
    % scatter((1:length(eVal)),sort(eVal, 'descend'), 10, "filled");
    % hold on;
    % yline([min_null max_null], '--', {'min' 'max'});
    % xlabel('Eigen value index');
    % ylabel('Eigen value');
    % box off;
    % axis padded;
    % ylim([-1.5 1])
    % % 
    % figure;
    % h = heatmap((cov_matrix), 'Colormap', parula);
    % h.YDisplayData = flipud(h.YDisplayData);
    % h.GridVisible = "off";
    % Labels = linspace(-stim_window*1e3,0,length(cov_matrix));
    % CustomLabels = string(Labels);
    % CustomLabels(mod(Labels,10) ~= 0) = " ";
    % h.XDisplayLabels = CustomLabels;
    % h.YDisplayLabels = flip(CustomLabels);
    % title("Spike Triggered Covariance (STC)");
    % 
    % figure;
    % h = heatmap((stim_prior_cov), 'Colormap', parula);
    % h.YDisplayData = flipud(h.YDisplayData);
    % h.GridVisible = "off";
    % Labels = linspace(-stim_window*1e3,0,length(stim_prior_cov));
    % CustomLabels = string(Labels);
    % CustomLabels(mod(Labels,10) ~= 0) = " ";
    % h.XDisplayLabels = CustomLabels;
    % h.YDisplayLabels = flip(CustomLabels);
    % title("Random stimulus prior");
    % 
    % figure;
    % h = heatmap(diff_cov, 'Colormap', parula);
    % h.YDisplayData = flipud(h.YDisplayData);
    % h.GridVisible = "off";
    % Labels = linspace(-stim_window*1e3,0,length(cov_matrix));
    % CustomLabels = string(Labels);
    % CustomLabels(mod(Labels,10) ~= 0) = " ";
    % h.XDisplayLabels = CustomLabels;
    % h.YDisplayLabels = flip(CustomLabels);
    % title("STC-C");
end