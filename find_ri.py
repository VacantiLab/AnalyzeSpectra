def find_ri(ic_smooth_dict,mz_vals,sat):
    #finds the retention indices

    import numpy as np
    import pdb

    sat_np = np.array(sat)

    #set define the alkane mz values, highest to lowest
    #    they range from 128 to 576, the last value in the series is not included below so it must be 114
    alkane_mz = np.arange(576,114,-14)
    alkane_nc = np.arange(40,7,-1) #the number of carbon atoms corresponding to the alkane mz's

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

    #for all of the recorded alkane mz ion count values, find the maximum
    max_alkane_val = np.amax(alkane_mz_maxv)

    #normalize the alkane mz-specific maximum ion count values by the overall alkane mz ion count maximum
    norm_mz_maxv = alkane_mz_maxv/max_alkane_val

    #create a matrix for better visualization during function creation
    alkane_matrix = np.matrix([alkane_mz_rec,norm_mz_maxv,alkane_mz_maxi])

    #defind a threshold for the normalized value above for considering an alkane mz present
    alk_ic_norm_thresh = 0.3

    #find the indices of those values above the threshold
    alkane_elut_ind = np.where(norm_mz_maxv > alk_ic_norm_thresh)[0]

    #convert the array containing the indices of sat vector for maximum readings at alkane mz's to intiger values
    alkane_mz_maxi = alkane_mz_maxi.astype(int)

    #find the number of carbons in the alkanes of those mz's above the threshold and convert to retention indices
    alkane_nc_rec = alkane_nc_rec[[alkane_elut_ind]]
    ri_rec = 100*alkane_nc_rec

    #find the indices corresponding to the scan number of those mz's above the threshold
    alkane_mz_maxi = alkane_mz_maxi[[alkane_elut_ind]]
    ri_sat = sat_np[[alkane_mz_maxi]]

    #Transform the sat array to an ri array

    #arrage retention index and corresponding retention time arrays in increasing order (should be in decreasing order previously)
    ri_rec = np.sort(ri_rec)
    ri_sat = np.sort(ri_sat)
    ri_array = np.array([]) #initialize retention index array

    for t in sat_np:
        test_array = ri_sat - t #find the difference of the retention times marking 100 ri units and the current retention time
        positive_inds = np.where(test_array >= 0)[0] #find all of the indices where the above differences are positive
        if len(positive_inds > 0): #if there are positive indices (you are not past the last time marking an increase of 100 ri units)
            k = positive_inds[0] #take the first positive index

            if k==0: #if you are below the first rt marking the rt of an alkane the interval used is the first interval
                ri_f = ri_rec[k+1]
                ri_i = ri_rec[k]
                rt_f = ri_sat[k+1]
                rt_i = ri_sat[k]

            if k > 0: #if you are above the first rt marking the rt of an alkane, the interval used flanks where you are
                ri_f = ri_rec[k]
                ri_i = ri_rec[k-1]
                rt_f = ri_sat[k]
                rt_i = ri_sat[k-1]

        #if there are no positive indices (you are past the last time marking an increase of 100 ri units)
        #    use the previous interval
        if len(positive_inds) == 0:
            last_ind = len(ri_rec) - 1
            ri_f = ri_rec[last_ind]
            ri_i = ri_rec[last_ind-1]
            rt_f = ri_sat[last_ind]
            rt_i = ri_sat[last_ind-1]

        #calculate the retention indices
        m = (ri_f - ri_i)/(rt_f-rt_i) #find the local slope
        r = t - rt_i #find the range
        ri_current = ri_i + m*r #calculate the current retention index
        ri_array = np.append(ri_array,ri_current) #update the array to store it in

    return(ri_array)
