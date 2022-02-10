def organize_ms_data(file_directory):
    import copy
    from pdb import set_trace#python debugger
    import importlib #allows fresh importing of modules
    from netCDF4 import Dataset #allows for opening netCDF files
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    from AnalyzeSpectra import remove_repeats #a custom function

    #read the netCDF file as a  netcdf4 object
        #ic = [ion counts for scan1 mz1, ion counts for scan1 mz2, ion counts for scan1 mz3, ....ion counts for scan107 mz1, ...]
        #mz = [mz for scan1 mz1, mz for scan1 mz2, mz for scan1 mz3, ...mz for scan107 mz1, ...]
        #si = [index where scan1 starts in ic and mz, index where scan2 starts in ic and mz, ...]
        #sat = [the scan acquisition time of scan1, the scan acquisition time of scan2, ...]
    ncdf = Dataset(file_directory,mode='r')
    ic = ncdf.variables['intensity_values'] #all of the ion count measurements for each scan. scans are concatenated into a single array
    mz = ncdf.variables['mass_values'] #all of the mz values for each of the 'intensity_values', all scans are concatenated as a single array
    si = ncdf.variables['scan_index'] #marks the python index of the starting position of each scan within the 'intensity_values'
    si = np.array(si)
    si = np.unique(si)
    sat = np.array(list(ncdf.variables['scan_acquisition_time'])) #the scan acquisition times corresponding to each scan (over mz values)
    sat_unique_values, sat_unique_indices = np.unique(sat, return_index=True)
        # somehow their can be two entries for the same scan acquistion time
        #     My guess is that the times are actually different, but rounding buy the machine makes them the same
        #     A probelm I ran into is that there would be more scan times than data sets for scan times, by a handful
    sat = sat_unique_values
    si = si[sat_unique_indices]
    tic = np.array(ncdf.variables['total_intensity']) #store the total ion count for each scan
    tic = tic[sat_unique_indices]
    n_scns = len(si) #the number of scans
    n_mz = len(mz) #the total number of recorded values
    mz_np = np.array(mz) #produces numpy float32 entries

    #create a data frame containing the mz values down the rows and the scan acquisition time values across the columns
    #   Older CDF files (from Christian's lab) contain extraneous scan information that needs to be removed.
    #   This is apparent because 'scan_index' has repeat values at the end
    if si[n_scns-1] == n_mz:
        n_scns = n_scns - 1
        sat = sat[0:n_scns]
        tic = tic[0:n_scns]

    #change mz values to type numpy float64 because that's what the library values are for mz and they need to agree for dictionary key  and pandas data frame row-name purposes
    mz_np2 = np.zeros(len(mz_np))
    for i in range(0,len(mz_np)):
        mz_np2[i] = np.float64(mz_np[i])
    mz_np = copy.copy(mz_np2)
    mz_np = np.around(mz_np,decimals=2) #round to have writable keys for dictionaries and pandas data frames

    mz_u = np.unique(mz_np) #unique scan values
    mz_u = np.sort(mz_u)

    #some measurements at mz's that round to the same whole mz appear in a single scan.
    #In these cases, the value closest to the whole mz value is retained, the other discarded
    #print('removing duplicate mz scans and updating scan indices for each SAT')
    #mz,ic,si = remove_repeats.remove_repeats(mz,ic,si)
    #The above is wrong - make a function here for binning based on the ms resolution

    #create a data frame containing the mz values down the rows and the scan acquisition time values across the columns
    ic_dct = dict()
    scan_indices = range(0,n_scns) #the scan indices range from 0 to 1 less thabn the # of scans because indexing starts at 0, so the 5th item has an index of 4
    for i in scan_indices:
        scan_start_index = si[i]
        if i < (n_scns-1): #if you are before the last scan index
            scan_end_index = si[i+1]
        if i == (n_scns-1): #if you at the last scan index
            scan_end_index = n_mz #this is not out of range because when specifiying a range to splice a list, the final value is not included

        # make a dictionary of pandas serieses
        #   each dictionary entry is labeled as the sat,
        #   contains the pandas series which is indexed by the mz values
        ics = ic[scan_start_index:scan_end_index]
        index_mzs = mz_np[scan_start_index:scan_end_index]

        # The mz values in a series can sometimes be duplicated - maybe due to instrument rounding
        #   remove the duplicate mzs and corresponding ic entries
        index_mzs_unique_values, index_mzs_unique_indices = np.unique(index_mzs, return_index=True)
        index_mzs = index_mzs_unique_values
        ics = ics[index_mzs_unique_indices]

        # Assemble the dictionary
        ic_dct[sat[i]] = pandas.Series(ics,index=index_mzs)


    #create the data frame with mz values down the rows and SATs across the columns from the dictionary of pandas serieses
    ic_df = pandas.DataFrame(ic_dct)

    #get an array of all mz values scanned (Agilent does not record all mz values for each scan because some mz values do not pass threshhold ioun counts in every scan)
    #mz_vals = ic_df.index.values

    #the number of mz values scanned is not necessarily the number of mz values within the scan range because some mz values may not have had any ion count values recorded.
    #I believe these values are not saved as 0 by the Agilent software, they are just discarded
    #n_mz = len(ic_df.iloc[:,0])

    #fill the NaN values with 0
    #NaN values exist because not all scans have measurements recorded for the same set of mz values (as explained above)
    ic_df = ic_df.fillna(value=0)

    return(ic_df,sat,n_scns,mz_u,tic)
