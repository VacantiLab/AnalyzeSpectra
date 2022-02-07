def GetFileBatch(file_directory,use_alkanes):
    # opens the batch_to_alkane.txt in the data location directory

    import pandas as pd
    import pdb
    import copy
    import numpy as np
    from os import listdir #needed to get a list of files
    import re #for regular expressions
    from pdb import set_trace

    #Get a list of all of the files in the specified directory
    files = listdir(file_directory)

    if 'file_to_batch.txt' in files:
        file_to_batch_map_directory = file_directory + 'file_to_batch.txt'
        file_name_to_batch = pd.read_table(file_to_batch_map_directory, sep="\t", header=0)
        batches = pd.unique(file_name_to_batch['batch'])

        batch_to_alkane_map_directory = file_directory + 'batch_to_alkane.txt'
        batch_to_alkane = pd.read_table(batch_to_alkane_map_directory, sep="\t", header=0)
        alkane_files = pd.unique(batch_to_alkane['alkane_file'])

    if 'file_to_batch.txt' not in files:
        #Remove filenames that are not NetCDF files
        netcdf_pattern = re.compile('.cdf$|.netcdf$',re.IGNORECASE)
            #creates a regex pattern matching '.cdf' or '.netcdf' at the end of a string (specified by $), where the case is not important
        filename_holder = copy.copy(files)
        for filename in filename_holder:
            is_netcdf = bool(re.search(netcdf_pattern,filename)) #determines if the pattern is in the filename
            if not is_netcdf:
                files.remove(filename)

        #place the files in alphabetical order so they are processed and printed in order
        files = sorted(files)

        #remove the alkane files from the file list
        if use_alkanes:
            files.remove('alkanes.CDF')

        #Make file to batch data frame
        file_name_to_batch_dic = {'file_name':files,'batch':['A']*len(files)}
        file_name_to_batch = pd.DataFrame(file_name_to_batch_dic)

        #batch to alkane data frame
        batch_to_alkane_dic = {'batch':['A'],'alkane_file':'alkanes.CDF'}
        batch_to_alkane = pd.DataFrame(batch_to_alkane_dic)

        #Make array of batches
        batches = np.array(['A'])

        #Make list of alkane files
        alkane_files = ['alkanes.CDF']

    return(file_name_to_batch,batch_to_alkane,alkane_files,batches)
