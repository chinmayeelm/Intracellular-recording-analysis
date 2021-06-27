function filtered_antennal_movement  = butter_filtfilt(data, fc, fs, order)
    
    %a= .9258; b=93.15; c=-1.455; %before 1 June 2021
    a=0.5334; b=516.5; c = -3.233; 
    antennal_movement = (b./(data -a)).^(1/3) + c;
%     [x, y] = butter(order, fc/(fs/2), 'low');
%     fvtool(x,y);

%     filtered_antennal_movement = filtfilt(x, y, antennal_movement);

    lpFilt = designfilt('lowpassiir','FilterOrder',order, ...
             'PassbandFrequency',fc,'PassbandRipple',0.05, ...
             'SampleRate',fs);
    filtered_antennal_movement = filtfilt(lpFilt, antennal_movement);
    
end
