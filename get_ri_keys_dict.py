def get_ri_keys_dict(ic_smooth_dict,ri_array,mz_vals):
    #inverts the dictionary containing mz values as the keys and ion counts for each retention index as the items
    #    a dictionary with retention indices as the keys and an array of ion counts for each mz is returned

    import numpy as np
    from pdb import set_trace

    #Initialize the dictionary holding the intensity vs. mz data for each timepoint
    ic_smooth_dict_timekeys = dict()
    k = 0
    for time_value in ri_array:
        ic_smooth_dict_timekeys[time_value] =  np.zeros(len(mz_vals))
    k = k+1

    j=0
    for plotted_mz in mz_vals:
        #fill in the time_key dictionary
        k = 0
        for time_value in ri_array:
            ic_smooth_dict_timekeys[time_value][j] =  ic_smooth_dict[plotted_mz][k]
            k = k+1
        j=j+1 #index of mz

    return(ic_smooth_dict_timekeys)
