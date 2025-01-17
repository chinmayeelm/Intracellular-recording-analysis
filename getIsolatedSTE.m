function STE = getIsolatedSTE(stimulus,raster, stim_window, fs)
%GETSTE This function returns all the spike triggering stimulus patterns
%   Stimulus pattern of stim_window length preceding every spike in the
%   raster is gathered.

rowNspikes = sum(raster,2);
raster_data = raster(rowNspikes~=0,:); 
[m,~] = size(raster_data);
STE = []; %nan([numSpikes stim_window* fs+1]);

parfor i=1:m
    spike_locs = find(raster_data(i,:)==1);
    valid_spike_locs = spike_locs(spike_locs >= stim_window*fs);
    isi = diff(valid_spike_locs);
    lowISIidx = find(isi<=4e-3*fs); %spikes closer than 4 ms.
    isolatedSpLocs = valid_spike_locs;
    isolatedSpLocs(lowISIidx+1) = [];
    spike_triggers = nan([length(isolatedSpLocs) stim_window*fs]);
    for j=1:length(isolatedSpLocs)
            
            spike_triggers(j,:) = stimulus(i,(isolatedSpLocs(j)-stim_window* fs+1):isolatedSpLocs(j));

        
    end
    STE = [STE; spike_triggers];

end
TF=isnan(STE(:,1));
STE(TF,:)=[];
end