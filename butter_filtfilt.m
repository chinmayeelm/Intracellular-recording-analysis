function filtered_antennal_movement  = butter_filtfilt(data, fc, fs, order,a, b, c)
    

    antennal_movement = (b./(data -a)).^(1/3) + c;
    [x, y] = butter(order, fc/(fs/2), 'low');
%     fvtool(x,y);

    filtered_antennal_movement = filtfilt(x, y, -antennal_movement); 
    % Multiplying -1 to the antennal movement because, 
    % "data", which is HES output voltage
    % is distance between the sensor and the magnet.
    % Above baseline voltage => Magnet closer to the sensor => Dorsad movement
    % Below baseline voltage => Magnet away from the sensor => Ventrad movement
    % But the calculation of antennal movement (distance) gives above baseline values
    % for movements away from sensor and below baseline values for
    % movements towards the sensor. The signal has to be flipped to
    % represent the movement correctly in terms of dorsad and ventrad.

    

%     lpFilt = designfilt('lowpassiir','FilterOrder',order, ...
%              'PassbandFrequency',fc,'PassbandRipple',0.05, ...
%              'SampleRate',fs);
%     filtered_antennal_movement = filtfilt(lpFilt, antennal_movement);
    
end
