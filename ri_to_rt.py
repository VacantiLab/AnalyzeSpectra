def ri_to_rt(sat_array,ri_array,ri):
    import numpy as np
    import pdb

    diff_array = ri_array-ri
    abs_array = np.absolute(diff_array)
    min_ind = np.argmin(abs_array)
    rt = sat_array[min_ind]
    return(rt)
