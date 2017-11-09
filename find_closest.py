def find_closest(value,array):
    import numpy as np
    import pdb

    dif_array = array - value
    dif_array_abs = np.absolute(dif_array)
    min_ind = np.argmin(dif_array_abs)
    min_val = array[min_ind]

    return(min_ind,min_val)
