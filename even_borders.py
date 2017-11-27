def even_borders(ic_smooth_dict,peak_start_i_dict,peak_end_i_dict,mz_vals):

    #move the lower peak border towards the center of the peak until it corresponds
    #    to an equal height of the other border
    #create a dictionary marking the peak indices where an index corresponds to a scan acquisition
    #the keys of the dictionary are mz values and the value is 0 for no peak at that index and 1 for a peak

    import numpy as np
    import copy
    import pdb

    #copy the input dictionaries so as not to change their values in the function
    peak_start_i_dict_copy = copy.deepcopy(peak_start_i_dict)
    peak_end_i_dict_copy = copy.deepcopy(peak_end_i_dict)

    #initialize the peak border dictionaries where the borders are of equal height
    peak_start_i_even_dict = dict()
    peak_end_i_even_dict = dict()

    #iterate through each mz value and relocate either the start or end border to make them equal height
    for mz in mz_vals:
        ion_counts = ic_smooth_dict[mz]
        peak_start_i_array = peak_start_i_dict_copy[mz]
        peak_end_i_array = peak_end_i_dict_copy[mz]
        n_peaks = len(peak_start_i_array)
        peak_iteration_array = np.arange(0,n_peaks)
        peak_start_i_even_array = copy.copy(peak_start_i_array)
        peak_end_i_even_array = copy.copy(peak_end_i_array)
        #iterate through each peak and determine if the start or end needs to be moved towards the center
        #    the one with a lower height will be moved inwards until the heights are equal
        for peak_iteration in peak_iteration_array:
            peak_start_less = False
            peak_start_i_even = peak_start_i_array[peak_iteration]
            peak_end_i_even = peak_end_i_array[peak_iteration]
            #if the beginning of the peak is lower, move it in
            while ion_counts[peak_start_i_even] < ion_counts[peak_end_i_even]:
                old_dif = ion_counts[peak_end_i_even] - ion_counts[peak_start_i_even]
                new_dif = ion_counts[peak_end_i_even] - ion_counts[peak_start_i_even + 1]
                #only accept the move if it brings the borders closer in height value
                accept_move = np.abs(new_dif) < np.abs(old_dif)
                if accept_move:
                    peak_start_i_even = peak_start_i_even + 1
                #if you can't get any closer by moving in the same direction, break the loop
                if not accept_move:
                    break
                #if you moved the start, indicate so so that the end isn't moved as well
                peak_start_less = True
            #if the end of the peak is lower, move it in
            while (ion_counts[peak_start_i_even] > ion_counts[peak_end_i_even]) & (not peak_start_less):
                old_dif = ion_counts[peak_end_i_even] - ion_counts[peak_start_i_even]
                new_dif = ion_counts[peak_end_i_even-1] - ion_counts[peak_start_i_even]
                #only accept the move if it brings the borders closer in height value
                accept_move = np.abs(new_dif) < np.abs(old_dif)
                if accept_move:
                    peak_end_i_even = peak_end_i_even - 1
                #if you can't get any closer by moving in the same direction, break the loop
                if not accept_move:
                    break
            peak_start_i_even_array[peak_iteration] = peak_start_i_even
            peak_end_i_even_array[peak_iteration] = peak_end_i_even
        peak_start_i_even_dict[mz] = peak_start_i_even_array
        peak_end_i_even_dict[mz] = peak_end_i_even_array

    #create the dictionary marking the peak indices where an index corresponds to a scan acquisition
    #    the keys of the dictionary are mz values and the value is 0 for no peak at that index and 1 for a peak

    #initialize the dictionary
    peak_range_dict = dict()
    n_scans = len(ion_counts)

    #create this dictionary one key value (mz) at a time
    for mz in mz_vals:
        print(mz)
        #initialize each key entry
        peak_range_dict[mz] = np.zeros(n_scans)
        all_indices = np.arange(0,n_scans)
        peak_start_i_even_array = peak_start_i_even_dict[mz]
        peak_end_i_even_array = peak_end_i_even_dict[mz]
        n_peaks = len(peak_start_i_even_array)
        peak_iterations = np.arange(0,n_peaks)
        #iterate through each scan to determine if it occurs during a peak elution
        for i in all_indices:
            #iterate through each peak and check to see if it is eluting during the current scan
            j = 0
            for p in peak_iterations:
                current_start = peak_start_i_even_array[p]
                current_end = peak_end_i_even_array[p]

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
