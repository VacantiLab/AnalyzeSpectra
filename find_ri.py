def find_ri(ic_smooth_dict,mz_vals):
    import numpy as np
    import pdb

    #set define the alkane mz values, highest to lowest
    #    they range from 128 to 576, the last value in the series is not included below so it must be 114
    alkane_mz = np.arange(576,114,-14)

    #for each mz vs time plot within the defined alkane mz values, find the maximum
    alkane_mz_rec = np.array([]) #initialize the recorded mz values array
    alkane_mz_maxv = np.array([]) #initialize the maximum value array
    alkane_mz_maxi = np.array([]) #initialize the maximum index array
    for i in alkane_mz:
        #check if the alkane mz has measured ion counts associated with it
        #    it could be outside the scan range or missed because its rt is too early or late
        if i in mz_vals:
            mz_v_t = ic_smooth_dict[i] #the mz vs. time plot for that mz
            max_val = np.amax(mz_v_t) #the maximum ion count for that mz
            max_ind = np.argmax(mz_v_t) #the index of the maximum ion count for that mz

            #update the initialized arrays
            alkane_mz_rec = np.append(alkane_mz_rec,i) #recorded mz's within the alkane mz values
            alkane_mz_maxv = np.append(alkane_mz_maxv,max_val) #the maximum ion count for that mz
            alkane_mz_maxi = np.append(alkane_mz_maxi,max_ind) #the index of that maximum value

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

    #find the mz values of those mz's above the threshold
    alkane_mz_rec = alkane_mz_rec[[alkane_elut_ind]]

    return(norm_mz_maxv)
