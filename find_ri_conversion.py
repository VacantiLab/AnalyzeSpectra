def find_ri_conversion(ic_smooth_dict,mz_vals,sat,coelut_dict_sat,coelut_dict_val_sat,sample_name):
    #finds the retention indices

    import numpy as np
    import pdb
    import importlib
    from AnalyzeSpectra import match_fingerprint
    from AnalyzeSpectra import get_alkane_dict
    importlib.reload(get_alkane_dict)

    sat_array = np.array(sat)

    #set define the alkane mz values, highest to lowest
    #    they range from 100 to 562, the last value in the series is not included below so it must be 114
    alkane_mz = np.arange(562,128,-14)
    alkane_nc = np.arange(40,9,-1) #the number of carbon atoms corresponding to the alkane mz's

    #for each mz vs time plot within the defined alkane mz values, find the maximum
    alkane_mz_rec = np.array([]) #initialize the recorded mz values array
    alkane_nc_rec = np.array([]) #initialize the recorded number of carbons array
    alkane_mz_maxv = np.array([]) #initialize the maximum value array
    alkane_mz_maxi = np.array([]) #initialize the maximum index array
    alkane_mz_sat = np.array([]) #initialize the sat array of the eluting alkane

    j = 0
    for i in alkane_mz:
        #check if the alkane mz has measured ion counts associated with it
        #    it could be outside the scan range or missed because its rt is too early or late
        if i in mz_vals:
            mz_v_t = ic_smooth_dict[i] #the mz vs. time plot for that mz
            max_val = np.amax(mz_v_t) #the maximum ion count for that mz
            max_ind = np.argmax(mz_v_t) #the index of the maximum ion count for that mz
            n_c = alkane_nc[j] #number of carbons in the current alkane

            #update the initialized arrays
            alkane_mz_rec = np.append(alkane_mz_rec,i) #recorded mz's within the alkane mz values
            alkane_nc_rec = np.append(alkane_nc_rec,n_c)
            alkane_mz_maxv = np.append(alkane_mz_maxv,max_val) #the maximum ion count for that mz
            alkane_mz_maxi = np.append(alkane_mz_maxi,max_ind) #the index of that maximum value
        j = j+1
    alkane_mz_maxi = alkane_mz_maxi.astype('int')

    #retrieve the dictionary containing the alkane peak profiles
    alkane_dict,alkane_names = get_alkane_dict.get_alkane_dict(alkane_nc_rec,alkane_mz_maxi,sat)

    #from here, use the alkane profile_dict to determine if the alkane is preent with fingerprint_match
    #    the instances of ri have already been replaced with rt in the necessary dictionaries
    #    the dictionaries available to this function should be passable to find_fingerprint

    #for each index where an alkane mz maximum was recorded, determine if that max is due to the alkane eluting
    #    if it is eluting, record the index where it elutes
    #    the recorded index is corresponds to the array containing all alkane masses that were scanned
    #        thus the sat, mz, and nc value of the alkanes corresponding can all be gathered from this array of indices
    #in progress
    alkane_elut_ind = np.array([])
    alkane_rt_array = np.array([])
    j = 0
    for alkane in alkane_names:
        metabolite = alkane
        ri_window = 1.5 #gives a little leway because peaks could be a little ugly if the chromatography is bad
        alkane_present,alkane_rt = match_fingerprint.match_fingerprint(sat,coelut_dict_sat,coelut_dict_val_sat,alkane_dict,mz_vals,ic_smooth_dict,metabolite,sample_name,ri_window)
        if alkane_present == True:
            alkane_elut_ind = np.append(alkane_elut_ind,j)
            alkane_rt_array = np.append(alkane_rt_array,alkane_rt)
        j = j+1
    alkane_elut_ind = alkane_elut_ind.astype('int')

    #find the number of carbons in the alkanes of those mz's above the threshold and convert to retention indices
    alkane_nc_rec = alkane_nc_rec[[alkane_elut_ind]]
    ri_rec = 100*alkane_nc_rec

    #find the indices corresponding to the scan number of those mz's above the threshold
    alkane_mz_maxi = alkane_mz_maxi[[alkane_elut_ind]]
    ri_sat = sat_array[[alkane_mz_maxi]]

    #Transform the sat array to an ri array

    #arrage retention index and corresponding retention time arrays in increasing order (should be in decreasing order previously)
    ri_rec = np.sort(ri_rec)
    ri_sat = np.sort(ri_sat)

    return(ri_sat,ri_rec)
