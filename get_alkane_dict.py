def get_alkane_dict(alkane_nc,alkane_mz_maxi,sat):
    import numpy as np
    import pdb

    alkane_dict = dict()

    alkane_names = np.array([])
    for nc in alkane_nc:
        nc = nc.astype('int')
        #to_append = np.array_str(nc) + '_carbon_alkane'
        to_append = nc.astype(str) + '_carbon_alkane'
        alkane_names = np.append(alkane_names,to_append )

    #upon EI fragmentation, alkanes have a peak at mz equal to their molecular weight and a tail peak at one above
    #    they also have another peak at mz = mw - 30 and a larger peak one above
    #the peak_profile dictionary created below represents this
    #    the magnitudes are provided to give the mw peak precedence for identification
    #    it does not mean this peak is larger (in fact it is not)

    j = 0
    for alkane in alkane_names:
        nc = alkane_nc[j]
        alkane_mw = 12*nc + 6 + 2*(nc-2)
        index_of_rt = alkane_mz_maxi[j].astype('int')
        rt_guess = sat[index_of_rt]
        alkane_dict[alkane] = dict()
        peak_profile = np.array([])
        peak_profile = np.append(peak_profile,alkane_mw)
        peak_profile = np.append(peak_profile,1e6)
        peak_profile = np.append(peak_profile,0.7)
        peak_profile = np.append(peak_profile,0.3)
        peak_profile = np.append(peak_profile,alkane_mw-30)
        peak_profile = np.append(peak_profile,9e5)
        peak_profile = np.append(peak_profile,0.3)
        peak_profile = np.append(peak_profile,0.7)
        alkane_dict[alkane]['peak_profile'] = peak_profile
        #the rt is given as an ri because peak matching on alkanes must be performed on rt
        #    they define the ri, so ri does not exist while searching for them
        #    ri is used to find the other peaks, thus this is called ri to hijack that function for finding alkanes
        alkane_dict[alkane]['ri'] = rt_guess
        j = j+1

    return(alkane_dict,alkane_names)
