raster_data = P(2).raster(1,:);
actual_stim = P(2).antennal_movement(1,:);
sta = P(2).STA(0.04*fs:end);
fs =10000;
t = linspace(0,15,length(raster_data));

locs = find(raster_data==1);
recon_stim_mat = zeros(length(locs),length(raster_data));
sta_len = length(sta);

for i = 1:length(locs)
    if (locs(i)-sta_len >= 0) 
        recon_stim_mat(i,(locs(i)-(sta_len-1)):locs(i)) = sta;
    end
end

recon_stim = sum(recon_stim_mat);
recon_stim = recon_stim-mean(recon_stim);
recon_stim_norm = (recon_stim-mean(recon_stim))/max(recon_stim);

actual_stim = actual_stim-mean(actual_stim);
actual_stim_norm = (actual_stim-mean(actual_stim))/max(actual_stim);
figure;
plot(t,actual_stim_norm); hold on;
plot(t,recon_stim_norm);
ylabel 'Antennal movement';
xlabel 'time (s)'
legend ('Actual stimulus','Reconstructed stimulus');

[r,lags] = xcorr(actual_stim,recon_stim);
figure();
stem(lags,r);