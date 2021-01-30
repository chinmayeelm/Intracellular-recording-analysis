
function [stim_order_sorted,stim_order_vector, idx] = sortfromnwb(nwb_in)
   
    idx = 0;
    stim = nwb_in.stimulus_presentation.get('mechanical_stimulus');
    stim_order = stim.stimulus_description;

    stim_order = string(stim_order);
    stim_order_vector = split(stim_order, ','); 

    [stim_order_sorted,idx] = sort(stim_order_vector);

end