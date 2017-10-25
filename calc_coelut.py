def calc_coelut(peak_ri_dict,mz_vals,ri_array):
    #calculates the coelution dictionary
    #    coelut_dict has each retention index as a key
    #    the corresponding items are arrays of the masses with peaks that elut at that ri
    #    the same peak can be considered eluting at neighboring ris because of the ri window

    import numpy as np

    #set the retention index window
    ri_window = 1.0

    #initialize the coelution dictionary
    coelut_dict = dict()

    #iterate through the retention indices
    for ri in ri_array:
        #initialize the array containing the mz values who have peaks at the current ri
        coelut_dict[ri] = np.array([])
        #iterate through the mz values to check for peaks
        for mz in mz_vals:
            #get a vector of the differences between peak ri's and the current ri at the current mz
            test_array = np.absolute(peak_ri_dict[mz] - ri)
            #determine if any of these ri's are within the window from the current ri
            test_logical = test_array < ri_window
            elutes = np.any(test_logical)
            #if the current mz elutes within an ri window of the current ri, append it to the elution vector for the current ri
            if elutes:
                coelut_dict[ri] = np.append(coelut_dict[ri],mz)

    return(coelut_dict)
