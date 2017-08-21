def ExtendBounds(vector,original_length):
    import pdb #python debugger
    import numpy as np #this is numpy, allows for data frame and matrix handling
    new_vec = np.zeros(5*len(vector))
    j=0
    for i in vector:
        new_vec[j] = i-2
        new_vec[j+1] = i-1
        new_vec[j+2] = i
        new_vec[j+3] = i+1
        new_vec[j+4] = i+2
        j=j+5

    too_small_indices = np.where(new_vec < 0)[0]
    too_large_indices = np.where(new_vec > (original_length-1))[0]
    to_remove = np.concatenate((too_small_indices,too_large_indices))
    if len(to_remove > 0):
        new_vec = np.delete(new_vec,to_remove)
    new_vec = np.unique(new_vec)
    new_vec = new_vec.astype(np.int64)

    return(new_vec)
