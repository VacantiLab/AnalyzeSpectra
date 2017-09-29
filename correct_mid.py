def correct_mid(mid_u,formula,CM_i):

    import importlib
    import pdb
    import numpy as np

    #the theoretical MIDs are the rows of the correction matrix.
    #    Their length corresponds to the number of rows of the right inverse of the correcion matrix
    n_theoretical_mid_entries = CM_i.shape[0]

    #the measured MID must have the same number of entries as each of the theoretical MIDs
    #    if it is short, add zeros to make up for the difference
    if len(mid_u) < n_theoretical_mid_entries:
        mid_u_appendage = np.zeros(n_theoretical_mid_entries-len(mid_u))
        mid_u = np.append(mid_u,mid_u_appendage)

    #calculate corrected MID
    mid_c = np.dot(mid_u,CM_i)
    mid_c = mid_c/sum(mid_c)

    return(mid_c)
