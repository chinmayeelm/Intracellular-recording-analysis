function [gcfrDecay, tDecay, phasicGCFR, tPhasic] = getPhasicDecayGCFR(P)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

[~,maxFRLoc] = max(P.avg_gcfr);

start_point = P.OFF_dur*P.fs+1;
stim_name_parts = split(P.stim_name);
delta_t = str2double(stim_name_parts(2));
rampEndIdx = start_point+delta_t*P.fs;

gcfrDecay = P.avg_gcfr(maxFRLoc:rampEndIdx);
phasicGCFR = P.avg_gcfr(start_point:rampEndIdx);


tDecay = linspace((1/P.fs), length(gcfrDecay)/P.fs, length(gcfrDecay));
tPhasic = linspace((1/P.fs), length(phasicGCFR)/P.fs, length(phasicGCFR));

end