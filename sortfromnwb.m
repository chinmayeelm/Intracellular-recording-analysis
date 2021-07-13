
function [stim_order_sorted,stim_order_vector, idx] = sortfromnwb(stim_order)
   
    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);

end