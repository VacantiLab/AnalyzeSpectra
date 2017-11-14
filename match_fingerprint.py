def match_fingerprint(ri_array,coelut_dict,coelut_dict_val):
    import numpy as np
    import importlib
    import pdb
    import find_closest
    importlib.reload(find_closest)

    #create the fingerprint dictionary
    #    each key is the first mz of a group of mz values
    #    the values are the fractional abundance of the first mz and subsequent integer increments
    fingerprint = dict()
    fingerprint[261] = np.array([0.059,0.016,0.022,0.698,0.137,0.061,0.076])
    fingerprint[233] = np.array([0.057,0.036,0.712,0.139,0.057])

    #define the retention index of the metabolite as found in the library
    metabolite_ri = 1388

    #define the window around the library retention index the metabolite is allowed to be found
    ri_window = 5
    ri_upper = metabolite_ri + ri_window
    ri_lower = metabolite_ri - ri_window

    #find the indices within the retention index array that border this retention index window
    metabolite_ri_ind,metabolite_ri = find_closest.find_closest(metabolite_ri,ri_array)
    ri_upper_ind,ri_upper = find_closest.find_closest(ri_upper,ri_array)
    ri_lower_ind,ri_lower = find_closest.find_closest(ri_lower,ri_array)

    #initialize a dictionary where the key value will be the mz groups (keys of fingerprint)
    #    the values for each key will be the retention indices where all mzs within that group are present
    metabolite_elut_ri_dict = dict()

    #get a list of the mz values heading each defined group of mz values
    group_mz_list = list(dict.keys(fingerprint))

    #initialize the values for each key of metabolite_elut_ri_dict to be an empty np array
    for mz in group_mz_list:
        metabolite_elut_ri_dict[mz] = np.array([])

    #iterate through the retention indices of the defined window
    #    at each retention index, the coeluting peaks are examined
    #    the set of coeluting peaks is tested against each group of mz values
    #    for each group of mz values, the retention indices where all members of the group of mz values are in the coeluting peaks is recoreded
    for i in range(ri_lower_ind,ri_upper_ind+1):
        ri = ri_array[i]
        peak_mz_array = coelut_dict[ri]
        peak_val_array = coelut_dict_val[ri]
        for mz in group_mz_list:
            n_mz = len(fingerprint[mz])
            mz_in_group = np.arange(mz,mz+n_mz)
            lib_in_samp_ind = np.isin(mz_in_group,peak_mz_array)
            lib_in_samp = np.all(lib_in_samp_ind)
            if lib_in_samp:
                metabolite_elut_ri_dict[mz] = np.append(metabolite_elut_ri_dict[mz],ri)
                #for each group of mz values, the key name is the first mz
                #the ri's where all members of the group are eluting is recorded

    #next objective is to find the times common to all members of metabolite_elut_ri_dict
    #also have a mechanism to say the metabolite is not present
    pdb.set_trace()
