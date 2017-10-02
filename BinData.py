def BinData(df,interval_val):
    #df: the data frame being binned
    #interval_val: the width of the bin interval

    import pdb
    import numpy as np
    import pandas
    import pdb

    #get the mz values scanned by the ms
    values = np.array(list(df.index.values))
    n_values = len(values)
    min_val = np.amin(values)

    #create an array containing the centers of the bins used to group the mz scanned values
    #    the interval_val cannot be smaller than the scan interval (that doesn't make sense anyways) because of
    #        how the length of the bin_centers array is set
    bin_centers = np.zeros(n_values)
    for i in range(0,n_values):
        bin_centers[i] = min_val + i*interval_val

    #map each actually scanned mz to a specified bin center
    for i in values:
        bin_diff = np.absolute(bin_centers-i)
        nearest_val_index = np.argmin(bin_diff)
        df.loc[i,'bin'] = np.round(bin_centers[nearest_val_index])

    #group the entries in 'bin' by summing the values with the same 'bin' sepcification
    df = df.groupby('bin').agg('sum')
    #    groups by the entry in bin and sums the values of the grouped rows

    return(df)
