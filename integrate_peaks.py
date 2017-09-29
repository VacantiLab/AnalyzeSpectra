def integrate_peaks(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,peak_start_i_dict,peak_end_i_dict,x_data_numpy,fragment_dict,fragments_list):
    import importlib #allows fresh importing of modules
    import pdb #python debugger
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    import scipy #contains simpsons integration function
    import fragment_library #a custom function
    importlib.reload(fragment_library) #reload the custom function in case it was changed
    import correct_mid
    importlib.reload(correct_mid)

    mzs_of_peak_starts = np.array(list(peak_start_t_dict.keys()))

    for frag_iter in fragments_list:
        mz = fragment_dict[frag_iter]['mz']
        rt = fragment_dict[frag_iter]['rt']

        for i in mz:
            #find the index of the peak in the mz curve (the first peak is 0, the second 1, the third 2, ...)
            possible_peak_starts = np.where(peak_start_t_dict[i] < rt)[0]
            if len(possible_peak_starts) > 0:
                prosp_peak_start_nm = max(possible_peak_starts)
            if len(possible_peak_starts) == 0:
                peak_present = False

            #find the time at which the peak is finished eluting
            if len(possible_peak_starts) > 0:
                prosp_peak_end_t = peak_end_t_dict[i][prosp_peak_start_nm]
                #if the peak ends after the retention time, then there is a peak (the the proposed peak was required to start before the retention time two lines ago)
                peak_present = rt < prosp_peak_end_t

            if not peak_present:
                fragment_dict[frag_iter]['areas'] = np.append(fragment_dict[frag_iter]['areas'],0)
            if peak_present:
                x_start_i = peak_start_i_dict[i][prosp_peak_start_nm]
                x_end_i = peak_end_i_dict[i][prosp_peak_start_nm]
                x_range_i = np.arange(x_start_i,x_end_i)
                x_range_t = x_data_numpy[x_range_i]
                y_range_ic = ic_smooth_dict[i][x_range_i]
                integrated_area = scipy.integrate.simps(y_range_ic,x_range_t)
                fragment_dict[frag_iter]['areas'] = np.append(fragment_dict[frag_iter]['areas'],integrated_area)
        fragment_dict[frag_iter]['tot_area'] = np.sum(fragment_dict[frag_iter]['areas'])
        if fragment_dict[frag_iter]['tot_area'] > 0:
            fragment_dict[frag_iter]['mid'] = fragment_dict[frag_iter]['areas']/fragment_dict[frag_iter]['tot_area']
        if fragment_dict[frag_iter]['tot_area'] == 0:
            fragment_dict[frag_iter]['mid'] = 0*fragment_dict[frag_iter]['areas']

        #correct the mids, will work if the MID is all zeros
        #    in order to print, all fragments must have a corrected mid
        #    these corrected MIDs must be the same length for a fragment across all samples (are formula dependent so they will be)
        mid_to_correct = fragment_dict[frag_iter]['mid']
        formula_of_mid_to_correct = fragment_dict[frag_iter]['formula']
        CM_i = fragment_dict[frag_iter]['CM_i']
        fragment_dict[frag_iter]['mid_c'] = correct_mid.correct_mid(mid_to_correct,formula_of_mid_to_correct,CM_i)

        #clever way to iterate through dictionary - not used here anymore
        #if fragment_dict[frag_iter]['tot_area'] == 0:
        #    fragment_dict[frag_iter]['mid'] = {k: 0 for k, v in fragment_dict[frag_iter]['areas'].items()}

    return(fragment_dict)
