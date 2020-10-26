def locate_overlap(ic_smooth_dict,peak_start_i_dict,peak_end_i_dict,mz_vals,peak_max_dict):

    #find the index borders designating the fractional-heights of each peak
    #    if the fractional height is designated as 0.5, the indices found designate the borders marking the half-height width
    #create a dictionary marking the peak overlap indices
    #    the keys of the dictionary are mz values and the value is 0 for no peak at that index and 1 for a peak

    import importlib #allows fresh importing of modules
    import pdb #python debugger
    import numpy as np
    import copy
    import pdb
    from AnalyzeSpectra import find_closest

    peak_height_fraction = 0.5

    #copy the input dictionaries so as not to change their values in the function
    peak_start_i_dict_copy = copy.deepcopy(peak_start_i_dict)
    peak_end_i_dict_copy = copy.deepcopy(peak_end_i_dict)

    #initialize the peak border dictionaries where the borders are of equal height
    peak_start_i_overlap_dict = dict()
    peak_end_i_overlap_dict = dict()

    #iterate through each mz value and relocate either the start or end border to make them equal height
    for mz in mz_vals:
        ion_counts = ic_smooth_dict[mz]
        peak_max_array = peak_max_dict[mz]
        peak_start_i_array = peak_start_i_dict_copy[mz]
        peak_end_i_array = peak_end_i_dict_copy[mz]
        n_peaks = len(peak_start_i_array)
        peak_iteration_array = np.arange(0,n_peaks)
        peak_start_i_overlap_array = copy.copy(peak_start_i_array)
        peak_end_i_overlap_array = copy.copy(peak_end_i_array)
        #iterate through each peak and determine if the start or end needs to be moved towards the center
        #    the one with a lower height will be moved inwards until the heights are equal
        for peak_iteration in peak_iteration_array:

            #find the maximum of the current peak
            peak_max = peak_max_array[peak_iteration]

            #find the border indices of the current peak
            peak_start_i = peak_start_i_array[peak_iteration]
            peak_end_i = peak_end_i_array[peak_iteration]

            #find the ion count values for the current peak
            peak_ic_array = ion_counts[peak_start_i:peak_end_i+1]

            #find the index of the max of the current peak
            peak_ic_array_max_i = np.argmax(peak_ic_array)

            #find the number of values in the current peak
            n_peak_points = len(peak_ic_array)

            #divide the peak into increasing and decreasing halves
            peak_ic_increasing = peak_ic_array[0:peak_ic_array_max_i]
            peak_ic_decreasing = peak_ic_array[peak_ic_array_max_i:n_peak_points]

            #find the ion counts that start and end the peak
            peak_start_ic = ion_counts[peak_start_i]
            peak_end_ic = ion_counts[peak_end_i]

            #calculate the height of the peak based on the distance of the max from the higher peak border
            height_bot_ref = max(peak_start_ic,peak_end_ic)

            #calculate the ion count that marks halfway to the max from the higher peak border
            half_height_ic = peak_height_fraction*(peak_max - height_bot_ref) + height_bot_ref

            #Initialize the indices which mark the region where the peak above the half_height_ic as defined above
            #    They must be initialized because the conditions described below may not be satisfied (due to artifacts) to search for them
            peak_start_i_overlap_array[peak_iteration] = peak_start_i
            peak_end_i_overlap_array[peak_iteration] = peak_end_i

            #Each of the increasing and decreasing vectors must contain values to search them
            #    Some may not because an artifact of baseline correcting is that few peaks have values of zero and could have length of 2 or maybe 1
            #    Baseline correcting is performed after peaks are found because it considers peak location
            #    Another artifact is that some "peaks" are always decreasing or increasing resulting in one of these vectors being empty
            #    Anything recorded as a peak that is actually an artifact could be removed (at an expense), but these values are usually in the "noise" regime
            if (len(peak_ic_increasing)>0) & (len(peak_ic_decreasing)>0):

                #Find the indices where the peak is closest to it's half_height_ic
                peak_ic_increasing_halfheight_i, peak_ic_increasing_halfheight_val = find_closest.find_closest(half_height_ic,peak_ic_increasing)
                peak_ic_decreasing_halfheight_i, peak_ic_decreasing_halfheight_val = find_closest.find_closest(half_height_ic,peak_ic_decreasing)

                #convert peak indices to overall indices for the mz value
                peak_start_halfheight_i = peak_start_i + peak_ic_increasing_halfheight_i
                peak_end_halfheight_i = peak_start_i + len(peak_ic_increasing) + peak_ic_decreasing_halfheight_i

                #store these half_height_ic interval demarcating indices for each peak in an array
                peak_start_i_overlap_array[peak_iteration] = peak_start_halfheight_i
                peak_end_i_overlap_array[peak_iteration] = peak_end_halfheight_i

        #store the arrays of half_height_ic interval demarcating indices for each mz in the dictionary initialized previously
        peak_start_i_overlap_dict[mz] = peak_start_i_overlap_array
        peak_end_i_overlap_dict[mz] = peak_end_i_overlap_array


    #create the dictionary marking the peak indices where an index corresponds to a scan acquisition
    #    the keys of the dictionary are mz values and the value is 0 for no peak at that index and 1 for a peak
    #    the values indicating a peak are assigned over a range
    #    say the peak has borders from indices 1000 to 1020
    #        the peak may be indicated to extend from the beginning of its half (or whatever the designated fraction is) height to the end of it's half-height

    #initialize the dictionary
    peak_range_dict = dict()
    n_scans = len(ion_counts)

    #create this dictionary one key value (mz) at a time
    for mz in mz_vals:
        #initialize each key entry
        peak_range_dict[mz] = np.zeros(n_scans)
        all_indices = np.arange(0,n_scans)
        peak_start_i_overlap_array = peak_start_i_overlap_dict[mz]
        peak_end_i_overlap_array = peak_end_i_overlap_dict[mz]
        n_peaks = len(peak_start_i_overlap_array)
        peak_iterations = np.arange(0,n_peaks)
        #iterate through each scan to determine if it occurs during a peak elution
        for i in all_indices:
            #iterate through each peak and check to see if it is eluting during the current scan
            j = 0
            for p in peak_iterations:
                current_start = peak_start_i_overlap_array[p]
                current_end = peak_end_i_overlap_array[p]

                #if the current scan falls before the current peak
                #    stop searching the peaks because subsequenc peaks have later starts (are in order)
                if i < current_start:
                    break

                #if the current scan is associated with a peak, mark it as a 1
                if (i >= current_start) & (i <= current_end):
                    peak_range_dict[mz][i] = 1
                    #once you find a peak associated with the current scan, you don't need to search the other peaks
                    break

                #if the current scan is past the current peak end
                #    that peak can be removed from consideration because the indices go in order
                if i > current_end:
                    new_peak_iterations_indices = np.arange(j+1,len(peak_iterations))
                    peak_iterations = peak_iterations[new_peak_iterations_indices]

                j = j+1

    return(peak_range_dict)
