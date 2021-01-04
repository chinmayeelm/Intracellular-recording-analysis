function [stim_order_vector, stim_order_sorted,idx] = sortfromtextfile(filename)
    txt_file = sprintf("%s.txt", filename);
    fileID = fopen(txt_file, 'r');
    stim_order = fscanf(fileID, '%s');

    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);
end