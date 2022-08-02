
function [stim_order_sorted, idx] = sortfromnwb(stim_order_vector)
   
%     stim_order = string(stim_order);
%     stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);

end