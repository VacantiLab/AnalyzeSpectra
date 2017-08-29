def BinData(df,bin_center):
    #df: the data frame being binned
    #bin_center: what you are centering on
    #    1 for centering on a whole number
    #    0.2 for centering to the nearest 0.2
    import numpy as np
    import pandas
    import pdb

    values = list(df.index.values)
    #n_values = len(values)
    #df.loc[:,'bin'] = values
    for i in values:
        df.loc[i,'bin'] = np.round(i/bin_center)*bin_center
        #    dividing by a decimal, rounding, and then multiplying by that decimal rounds to that nearest decimal

    df = df.groupby('bin').agg('sum')
    #    groups by the entry in bin and sums the values of the grouped rows

    return(df)
