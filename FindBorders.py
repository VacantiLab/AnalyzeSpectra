def FindBorders(indexes,y_base_cor,sat):
    import pdb
    import numpy as np

    borders_prior_indices = np.zeros(len(indexes))
    borders_after_indices = np.zeros(len(indexes))

    #for each peak index, find the border prior to and following the corresponding peak
    j=0
    for i in indexes:
        #for each peak, the initial border candidates are 50% of the peak height
        border_candidates = 0.5*y_base_cor[i] >= y_base_cor

        #find the initial beginning border by finding the nearest border candidate prior to the peak index
        bor_cand_prior = border_candidates[0:i] #find all of the border candidates prior to the peak
        bor_cand_prior[0] = True #the first value is true because if all else fails the peak must begin where the data begins
        bor_can_prior_index = np.where(bor_cand_prior) #get the indices of the True values
        border_prior_index = np.max(bor_can_prior_index) #find the maximum index of True values prior to the peak

        #keep searching before until you increase in value
        values_prior = y_base_cor[0:(border_prior_index+1)] #get the values prior to the border index that was just found
        to_subtract = 0 #initialize the value you will subtract from that index
        #if the previous value is less than the value at the beginning border index, the index of that value becomes the new opening border index
        for k in range(1,len(values_prior)-1):
            if values_prior[-k]<=values_prior[-(k+1)]: # -k indexes the kth value from the end of the string
                break
            to_subtract = to_subtract+1
        borders_prior_indices[j] = border_prior_index - to_subtract #update the beginning border index

        #find the initial ending border by finding the nearest border candidate after the peak index
        bor_can_after = border_candidates[(i+1):len(border_candidates)] #find all of the border candidates following to the peak
        bor_can_after[len(bor_can_after)-1] = True #the last value is true because if all else fails the peak must end where the data ends
        bor_can_after_index = np.where(bor_can_after) #get the indices of the True values
        border_after_index = i + np.min(bor_can_after_index) #find the minimum index of True values prior to the peak

        #keep searching after until you increase in value
        values_after = y_base_cor[border_after_index:len(y_base_cor)] #get the values following the border index that was just found
        to_add = 0 #initialize the value you will add to that index
        #if the next value is less than the value at the end border index, the index of that next value becomes the new end border index
        for k in range(0,len(values_after)-1):
            if values_after[k]<=values_after[k+1]:
                break
            to_add = to_add+1
        borders_after_indices[j] = border_after_index + to_add #update the end border index

        #if the beginning of a peak is before the location of the previous peak
            #the beginning of that peak is moved to the minimum location after the previous peak and before the current peak
        if j > 0:
            current_peak_index = indexes[j]
            prior_peak_index = indexes[j-1]
            prior_border_index = borders_prior_indices[j]
            if prior_border_index < prior_peak_index:
                local_min_index = np.argmin(y_base_cor[prior_peak_index:current_peak_index+1])
                borders_prior_indices[j] = prior_peak_index + local_min_index

        #if the end of a peak is after the location of the next peak
            #the end of that peak is moved to the minimum location between the current and next peak
        if j < len(indexes)-1:
            current_peak_index = indexes[j]
            next_peak_index = indexes[j+1]
            after_border_index = borders_after_indices[j]
            if next_peak_index < after_border_index:
                local_min_index = np.argmin(y_base_cor[current_peak_index:next_peak_index+1])
                borders_after_indices[j] = current_peak_index + local_min_index

        j = j+1
    borders_prior_indices = borders_prior_indices.astype(int)
    borders_after_indices = borders_after_indices.astype(int)
    return(borders_prior_indices,borders_after_indices)
