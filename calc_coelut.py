def calc_coelut(peak_ri_dict,mz_vals,ri_array,ic_smooth_dict,peak_overlap_dictionary):
    #calculates the coelution dictionary
    #    coelut_dict has each retention index as a key
    #    the corresponding items are arrays of the masses with peaks that elut at that ri
    #    the same peak can be considered eluting at neighboring ris because of the ri window

    import numpy as np

    #set the retention index window
    ri_window = 1.0

    #initialize the coelution dictionary and the dictionary containing the values of corresponding peaks
    #    the keys of the two dictionaries are retention indices
    coelut_dict = dict()
    coelut_dict_val = dict()

    #iterate through the retention indices
    ri_i = 0
    for ri in ri_array:
        #initialize the array containing the mz values who have peaks at the current ri, and the array of values at those peaks
        coelut_dict[ri] = np.array([])
        coelut_dict_val[ri] = np.array([])
        #iterate through the mz values to check for peaks
        for mz in mz_vals:
            elutes = False
            if peak_overlap_dictionary[mz][ri_i] == 1:
                elutes = True

                ##THE OLD WAY GIVEN AN RI WINDOW
                ##get a vector of the differences between peak ri's and the current ri at the current mz
                #test_array = np.absolute(peak_ri_dict[mz] - ri)
                ##determine if any of these ri's are within the window from the current ri
                #test_logical = test_array < ri_window
                #elutes = np.any(test_logical)


            #if the current mz elutes within an ri window of the current ri, append it to the elution vector for the current ri
            #    also append the value of that peak at the current ri to the respective array
            if elutes:
                coelut_dict[ri] = np.append(coelut_dict[ri],mz)
                peak_value = ic_smooth_dict[mz][ri_i]
                coelut_dict_val[ri] = np.append(coelut_dict_val[ri],peak_value)
        ri_i = ri_i + 1

    return(coelut_dict,coelut_dict_val)
