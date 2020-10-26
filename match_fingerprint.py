def match_fingerprint(ri_array,coelut_dict,coelut_dict_val,metabolite_dict,mz_vals,ic_smooth_dict,metabolite,sample_name,ri_window):
    import numpy as np
    import importlib
    import pdb
    from AnalyzeSpectra import find_closest

    #find the mz index range of the scan
    mz_scan_start = mz_vals[0]
    mz_scan_end = mz_vals[len(mz_vals)-1]

    #define the retention index of the metabolite as found in the library
    metabolite_ri = metabolite_dict[metabolite]['ri']

    #determine the metabolite fingerprint array
    metabolite_peak_profile = metabolite_dict[metabolite]['peak_profile']

    #convert the peak profile stored in the library into the metabolite fingerprint and group total ion count dictionaries used by this function
    fingerprint = dict()
    group_tic_dict = dict()
    j = 0
    just_passed_key = False
    for item in metabolite_peak_profile:
        if (item > 1) & (not just_passed_key):
            current_key = item
            just_passed_key = True
            fingerprint[current_key] = np.array([])
            group_tic_dict[current_key] = metabolite_peak_profile[j+1]
        if item <= 1:
            fingerprint[current_key] = np.append(fingerprint[current_key],item)
            just_passed_key = False
        j = j+1

    #test a fingerprint
    #fingerprint = dict()
    #fingerprint[148] = np.array([0.3,0.4,0.2,0.1])
    #fingerprint[233] = np.array([0.2,0.4,0.3,0.1])
    #fingerprint[647] = np.array([0.2,0.4,0.4])
    #fingerprint[645] = np.array([0.4,0.5,0.05,0.03,0.02])
    #fingerprint[648] = np.array([0.1,0.4,0.45,0.05])
    #fingerprint[651] = np.array([0.3,0.5,0.15,0.05])
    #fingerprint[345] = np.array([0.2,0.4,0.3,0.1])
    #fingerprint[445] = np.array([0.2,0.4,0.3,0.1])
    #fingerprint[456] = np.array([0.2,0.4,0.3,0.1])

    #test a corresponding group_tic_dict
    #group_tic_dict = dict()
    #group_tic_dict[148] = 1409875
    #group_tic_dict[233] = 7890523
    #group_tic_dict[647] = 290856
    #group_tic_dict[645] = 8943679
    #group_tic_dict[648] = 7812520
    #group_tic_dict[651] = 10098763
    #group_tic_dict[345] = 134985
    #group_tic_dict[445] = 495863254
    #group_tic_dict[456] = 36846

    #trim the peak profile so that entries outiside of the mz scan range are not considered in the match
    #    if an mz key is below the mz scan range, remove member entries that are below the scan range and rename the key of the group
    #    if the mz key is above the mz scan range, remove it and all members
    #    remove all members of any mz group which correspond to mz values above the scan range
    fingerprint,group_tic_dict = trim_peak_profile(fingerprint,group_tic_dict,mz_scan_start,mz_scan_end)

    #rank the groups by abundance and only retain those above a specified percentile
    #    this allows for the lower abundance dictionary fragments to be missing in the spectrum and still have a match
    #    the percentile excluded depends on the number of groups there are in a library
    ranked_fingerprint,ranked_group_tic_dict = return_ranked_fingerprint(fingerprint,group_tic_dict)
    #define the window around the library retention index the metabolite is allowed to be found
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
    group_mz_list = list(dict.keys(ranked_fingerprint))

    #get a list of all group mz values: not just the ones leading and naming the group
    group_mz_vals_all = np.array([])
    for mz_group_name in group_mz_list:
        group_mz_vals_all = np.append(group_mz_vals_all,mz_group_name)
        mz_group_length = len(ranked_fingerprint[mz_group_name])
        to_add_then_append_array = np.arange(1,mz_group_length+1)
        for i in to_add_then_append_array:
            to_append = mz_group_name + i
            group_mz_vals_all = np.append(group_mz_vals_all,to_append)

    #initialize the values for each key of metabolite_elut_ri_dict to be an empty np array
    for mz in group_mz_list:
        metabolite_elut_ri_dict[mz] = np.array([])

    #iterate through the retention indices of the defined window
    #    at each retention index, the coeluting peaks are examined
    #    the set of coeluting peaks is tested against each group of mz values
    #    for each group of mz values, the retention indices where all members of the group of mz values are in the coeluting peaks is recoreded
    ri_index_range = range(ri_lower_ind,ri_upper_ind+1)
    for i in ri_index_range:
        ri = ri_array[i]
        peak_mz_array = coelut_dict[ri]
        peak_val_array = coelut_dict_val[ri]
        #see if the groups specified by the leading mz's are considered eluting at the current ri
        #    store these in a dictionary where the mz values specifying the groups are the keys
        #        the values are the retention indices where the group is considered present
        for mz in group_mz_list:
            group_present = check_group(mz,ranked_fingerprint,peak_mz_array,peak_val_array,ri,mz_scan_end,metabolite,sample_name)
            if group_present:
                metabolite_elut_ri_dict[mz] = np.append(metabolite_elut_ri_dict[mz],ri)
                #for each group of mz values, the key name is the first mz
                #the ri's where group-eluting criteria is satisfied is recorded

    #find the retention indices where all groups (as filtered above) are present
    j = 0
    for mz in group_mz_list:
        if j==0:
            intersection = metabolite_elut_ri_dict[mz]
        if not j==0:
            intersection = np.intersect1d(intersection,metabolite_elut_ri_dict[mz])
        j = j+1

    #find the ris where the mzs of the groups are at a max
    #    the median of these is reported as the retention index of the metabolite
    metabolite_present = False #initialize this value
    metabolite_retention_index = np.array([]) #initialize this value
    if len(intersection) >= 1: #there must be 1 point of overlap
        max_ri_array = np.array([])
        max_value_array = np.array([])
        for mz_to_max in group_mz_vals_all:
            intersection_indices = np.array([])
            for i in intersection:
                intersection_index = np.where(ri_array==i)
                intersection_indices = np.append(intersection_indices,intersection_index)
                intersection_indices = intersection_indices.astype(int)
            value_array = ic_smooth_dict[mz_to_max][intersection_indices]
            max_index = np.argmax(value_array)
            max_ri_index = intersection_indices[max_index]
            max_ri = ri_array[max_ri_index]
            max_ri_array = np.append(max_ri_array,max_ri)

            max_value = np.amax(value_array)
            max_value_array = np.append(max_value_array,max_value)

        max_ri = np.median(max_ri_array)
        max_ri_ind,max_ri = find_closest.find_closest(max_ri,ri_array) #taking a median one line abouve may make the ri a value outside of the ri_array. it must be in the ri_array for indexing purposes.

        max_value = np.amax(max_value_array)

        if sample_name != 'alkanes':
            metabolite_present = True
            metabolite_retention_index = max_ri

        # if i am looking at the alkanes sample, the first mz value must have a max intensity meeting a certain threshold for it to be considered present
        # this tends to exclude larger alkanes that did not actually elute but have negligible peaks somewhere (I don't know why)
        max_value_threshold = 500
        alkane_parent_ion = group_mz_vals_all[0]
        alkane_parent_ion_max_val = max_value_array[0]
        if sample_name == 'alkanes':
            threshold_condition = alkane_parent_ion_max_val > max_value_threshold
            if threshold_condition:
                metabolite_present = True
                metabolite_retention_index = max_ri

            # the intensity of the next alkane parent ion must also not be above the threshold if it is present
            # this excludes finding a smaller alkane that is not there because it is a fragment of a larger alkane
            next_alkane_mz = alkane_parent_ion + 14
            co_eluting_mz_array = coelut_dict[metabolite_ri] # if you use max_ri this is a problem for some reason
            co_eluting_mz_val_array = coelut_dict_val[metabolite_ri]
            if np.isin(next_alkane_mz,co_eluting_mz_array):
                next_alkane_i = np.where(next_alkane_mz == co_eluting_mz_array)[0]
                if co_eluting_mz_val_array[next_alkane_i] > max_value_threshold:
                    metabolite_present = False

    return(metabolite_present,metabolite_retention_index)

#Supporting Functions
################################################################################
################################################################################
def trim_peak_profile(fingerprint,group_tic_dict,mz_scan_start,mz_scan_end):

    import numpy as np
    import pdb

    #retrieve a list of all mz values beginning a group
    group_mz_list = list(dict.keys(fingerprint))

    #iterate througn those values
    for mz in group_mz_list:
        #keep track if the mz key is removed or renamed so a removed key isn't referenced later
        mz_key_removed_renamed = False
        #if you find an mz value beginning a group that is smaller than the smallest mz value scanned
        #    remove the key and all associated values if all associated values are less than the first scanned mz value
        #    if the mz key value (group start value) is less than the scan start but there are mz values within the group that are not
        #        remove all values smaller than the mz scan start, renormalize the group, and rename the group (the key) after the new smallest value
        if mz < mz_scan_start:
            group_mz = fingerprint[mz] #retrieve an array of all of the mz values in that group
            n_group_mz = len(group_mz) #retrieve the number of mz values in that group
            group_tic = group_tic_dict[mz] #retrieve the total ion count associated with that group
            #iterate through each individual value in that group
            #keep track of the number of values that must be removed (those less than the first mz value scanned)
            remove_count = 0
            for i in np.arange(0,n_group_mz):
                if mz+i < mz_scan_start:
                    remove_count = remove_count + 1
            #remove the values corresponding to thouse found to be smaller than the smallest value scnaned
            indices_to_delete = np.arange(0,remove_count)
            tic_to_subtract = np.sum(group_mz[indices_to_delete]*group_tic)
            group_mz = np.delete(group_mz,indices_to_delete)
            group_tic = group_tic - tic_to_subtract
            #update the dicionary key entry appropriately
            if len(group_mz) == 0: #remove the key entirely if it corresponds to an empty array
                del fingerprint[mz]
                del group_tic_dict[mz]
                mz_key_removed_renamed = True
            if len(group_mz) > 0: #rename the entry after the smallest (first) value
                fingerprint[mz] = group_mz/np.sum(group_mz)
                group_tic_dict[mz] = group_tic
                new_key = mz + remove_count
                fingerprint[new_key] = fingerprint.pop(mz)
                group_tic_dict[new_key] = group_tic_dict.pop(mz)
                mz_key_removed_renamed = True

        #if the current mz key is higher than the highest mz scanned, remove it from the metabolite fingerprint
        if (mz > mz_scan_end) & (not mz_key_removed_renamed):
            del fingerprint[mz]
            del group_tic_dict[mz]
            mz_key_removed_renamed = True

        #if the current mz key is lower than the highest mz scanned
        #    check to see if any of the members of that group are higher than the highest mz scanned
        #    remove them from the metabolite profile if they are
        if (mz <= mz_scan_end) & (not mz_key_removed_renamed):
            group_mz = fingerprint[mz]
            group_tic = group_tic_dict[mz] #retrieve the total ion count associated with that group
            n_group_mz = len(group_mz)
            indices_to_delete = np.array([])
            j = 0
            for i in np.arange(0,n_group_mz):
                if mz + i > mz_scan_end:
                    indices_to_delete = np.append(indices_to_delete,j)
                j = j+1
            indices_to_delete = indices_to_delete.astype(int)

            if len(indices_to_delete > 0):
                tic_to_subtract = np.sum(group_mz[indices_to_delete]*group_tic)
                group_mz = np.delete(group_mz,indices_to_delete)
                group_tic = group_tic - tic_to_subtract
                fingerprint[mz] = group_mz/np.sum(group_mz)
                group_tic_dict[mz] = group_tic

    return(fingerprint,group_tic_dict)

################################################################################
################################################################################

def return_ranked_fingerprint(fingerprint,group_tic_dict):
    import numpy as np
    import pdb
    import copy

    #copy the fingerprint and group_tic_dict so as not to change their values in this function
    ranked_fingerprint = copy.copy(fingerprint)
    ranked_group_tic_dict = copy.copy(group_tic_dict)

    #retrieve a list of all mz values beginning a group
    group_mz_list = list(dict.keys(fingerprint))


    #set the percentile above which the groups are retained (based on abundance)
    n_groups = len(group_mz_list)
    if n_groups <= 3:
        percentile = 0
    if (n_groups > 3) & (n_groups <= 5):
        percentile = 20
    if (n_groups > 5) & (n_groups <= 10):
        percentile = 40
    if (n_groups > 10):
        percentile = 50

    if percentile > 0:
        #retrieve a list of corresponding total ion counts
        ion_counts = np.array([])
        for mz in group_mz_list:
            ion_counts = np.append(ion_counts,group_tic_dict[mz])

        #retrieve the abundance marking the cut-off percentile
        cut_off_abundance = np.percentile(ion_counts,percentile)

        #assemble a list of the mz values marking groups not to consider
        mzs_to_remove = np.array([])
        for mz in group_mz_list:
            if group_tic_dict[mz] < cut_off_abundance:
                mzs_to_remove = np.append(mzs_to_remove,mz)

        #remove the above assembled list of mz values from the considered library values
        for mz in mzs_to_remove:
            del ranked_fingerprint[mz]
            del ranked_group_tic_dict[mz]

    #return the updated values
    return(ranked_fingerprint,ranked_group_tic_dict)

################################################################################
################################################################################

def check_group(mz,ranked_fingerprint,peak_mz_array,peak_val_array,ri,mz_scan_end,metabolite,sample_name):
    #this function sets the rules to determine if a group led by an mz is present
    import numpy as np
    import pdb

    #Initialize the three specifcations that determine if the group is present
    quantity_requirement = False
    quantity_requirement_with_shift = False
    magnitude_requirement = False

    #specify the mz values contained within the group
    n_mz = len(ranked_fingerprint[mz])
    mz_in_group = np.arange(mz,mz+n_mz)

    #for the mz values in this group, determine how many of them are in the coelution array
    #    the coelution array is supplied for a given retention index
    #if enough of the mz values are coeluting, the quantity requirement is met
    index_in_coelut_array = np.array([])
    for i in mz_in_group:
        index_in_coelut = np.where(peak_mz_array==i)
        index_in_coelut_array = np.append(index_in_coelut_array,index_in_coelut)
    if len(index_in_coelut_array) >= n_mz-1:
        quantity_requirement = True

    #metabolites may labeled, so the mz group is extended and a new quantity requirement is set
    #    this new quantity requirement OR the previous one must be satisfied to say the metabolite is present
    index_in_coelut_array = np.array([])
    last_mz_in_group = mz_in_group[n_mz-1]
    mzs_to_append = np.array([last_mz_in_group+1,last_mz_in_group+2,last_mz_in_group+3])
    mz_group_extended = np.append(mz_in_group,mzs_to_append)
    n_mz_extended = len(mz_group_extended)
    for i in mz_group_extended:
        index_in_coelut = np.where(peak_mz_array==i)
        index_in_coelut_array = np.append(index_in_coelut_array,index_in_coelut)
    if len(index_in_coelut_array) >= n_mz_extended-2:
        quantity_requirement_with_shift = True

    #the value associated with the mz prior to that defining the group cannot be larger than
    #    all of the first 4 values of the group
    #    you must also consider if you are at the end of the mz scan range
    #        if so this test cannot be used
    if not mz-1 in peak_mz_array:
        magnitude_requirement = True
    at_end_mz_scan_range = mz+4 > mz_scan_end

    if (mz-1 in peak_mz_array) & (not at_end_mz_scan_range):
        prior_index = np.where(peak_mz_array==mz-1)[0]
        prior_value = peak_val_array[prior_index]
        group_indices = np.arange(prior_index+1,prior_index+5)
        n_in_group_val_array = 0
        prior_value_smaller_array = np.array([])
        n_prior_value_smaller = 0
        for i in group_indices:
            #at the given ri, all peaks may not be eluting
            if len(peak_val_array) >= i+1:
                group_value = peak_val_array[i]
                n_in_group_val_array = n_in_group_val_array + 1
                prior_value_smaller = prior_value < group_value
                prior_value_smaller_array = np.append(prior_value_smaller_array,prior_value_smaller)
        n_prior_value_smaller = sum(prior_value_smaller_array)

        #if there is a single value the prior magnitude is less than in the group
        #    then this criteria is met
        if (n_prior_value_smaller > 0):
            magnitude_requirement = True

    #if you are at the end of the scan range, this criteria is not applicable and is given to be satisfied
    if at_end_mz_scan_range:
        magnitude_requirement = True

    #determine if the group is present from the above three specified criteria
    group_present = (quantity_requirement | quantity_requirement_with_shift) & magnitude_requirement

    return(group_present)
