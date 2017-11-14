def convert_rt_ri(ri_sat,ri_rec,sat_array):
    import numpy as np
    import pdb

    #initialize retention index array
    ri_array = np.array([])

    #iterate through each scan acquisition time and convert it to a retention index
    for t in sat_array:
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
