def integrate_peaks(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,peak_start_i_dict,peak_end_i_dict,x_data_numpy):
    import importlib #allows fresh importing of modules
    import pdb #python debugger
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    import scipy #contains simpsons integration function
    import fragment_library #a custom function
    importlib.reload(fragment_library) #reload the custom function in case it was changed

    fragment_dict = fragment_library.fragment_library()
    fragments_ls = list(fragment_dict.keys())

    for frag_iter in fragments_ls:
        mz = fragment_dict[frag_iter]['mz']
        rt = fragment_dict[frag_iter]['rt']

        for i in mz:
            #find the index of the peak in the mz curve (the first peak is 0, the second 1, the third 2, ...)
            prosp_peak_start_nm = max(np.where(peak_start_t_dict[i] < rt)[0])
            #find the time at which the peak is finished eluting
            prosp_peak_end_t = peak_end_t_dict[i][prosp_peak_start_nm]
            #if the peak ends after the retention time, then there is a peak (the the proposed peak was required to start before the retention time two lines ago)
            peak_present = rt < prosp_peak_end_t
            if not peak_present:
                fragment_dict[frag_iter]['areas'][i] = 0.0
            if peak_present:
                x_start_i = peak_start_i_dict[i][prosp_peak_start_nm]
                x_end_i = peak_end_i_dict[i][prosp_peak_start_nm]
                x_range_i = np.arange(x_start_i,x_end_i)
                x_range_t = x_data_numpy[x_range_i]
                y_range_ic = ic_smooth_dict[i][x_range_i]
                fragment_dict[frag_iter]['areas'][i] = scipy.integrate.simps(y_range_ic,x_range_t)
        fragment_dict[frag_iter]['tot_area'] = sum(fragment_dict[frag_iter]['areas'].values())
        fragment_dict[frag_iter]['mid'] = {k: v / fragment_dict[frag_iter]['tot_area'] for k, v in fragment_dict[frag_iter]['areas'].items()}

    return(fragment_dict)
