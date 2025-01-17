function STE = getSTE(stimulus,raster, stim_window, fs)
%GETSTE This function returns all the spike triggering stimulus patterns
%   Stimulus pattern of stim_window length preceding every spike in the
%   raster is gathered.

rowNspikes = sum(raster,2);
raster_data = raster(rowNspikes~=0,:); 
[m,~] = size(raster_data);
STE = []; %nan([numSpikes stim_window* fs+1]);

for i=1:m
    spike_locs = find(raster_data(i,:)==1);
    valid_spike_locs = spike_locs(spike_locs >= stim_window*fs);
    spike_triggers = nan([length(valid_spike_locs) stim_window*fs]);
    for j=1:length(valid_spike_locs)
            
            spike_triggers(j,:) = stimulus(i,(valid_spike_locs(j)-stim_window* fs+1):valid_spike_locs(j));

        
    end
    STE = [STE; spike_triggers];

end
TF=isnan(STE(:,1));
STE(TF,:)=[];
end