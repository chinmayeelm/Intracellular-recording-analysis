function crossing = schmitt_trigger(signal, low_level, high_level)
%SCHMITT_TRIGGER 
%   
find(diff(signal > low_level & signal < high_level))
end

