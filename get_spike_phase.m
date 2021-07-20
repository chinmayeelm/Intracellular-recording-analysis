function spike_phase_val = get_spike_phase(index, fs, period)

spike_phase_val = (index/fs)*((2*pi)/period);

end