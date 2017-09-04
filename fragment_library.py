def fragment_library():

    import pdb
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas
    import re

    file_name_read = '/Users/Nate/Desktop/netcdf_test/tbdms_lib.txt'
    fragment_dict = dict()

    #open a .txt file with the fragment information and import it
    #.txt file has format:
    #    fragment: pyr174
    #    formula: C6H12N1O3Si1
    #    rt: 450
    #    mz: 174 175 176 177 178
    #
    #    fragment: lac233
    #    ...
    with open(file_name_read, 'r') as read_file:
        #read through the lines in the file one by one
        outer_read_line = 0
        for line in read_file:
            outer_read_line = outer_read_line + 1 #iterate the line number
            line_split = line.split(':') #split the line into a list of strings, colons mark separations
            #if the first word on the line is 'fragment', gather the fragment information
            if line_split[0]=='fragment':
                fragment_name = line_split[1].lstrip().rstrip() #remove white space characters from the left and right of the fragment name
                fragment_dict[fragment_name] = dict() #initialize a dictionary for the current fragment
                fragment_dict[fragment_name]['areas'] = dict() #initialize the peak areas dictionary
                fragment_dict[fragment_name]['mid'] = dict() #initialize the MID dictionary
                #read through the same file line by line, but starting from the beginning
                inner_read_line = 0
                with open(file_name_read, 'r') as inner_read_file:
                    for inner_line in inner_read_file:
                        inner_read_line = inner_read_line + 1 #iterate the line number
                        inner_line_split = inner_line.split(':') #split the line into a list of strings
                        #read the fragment formula
                        if inner_read_line == outer_read_line + 1:
                            fragment_formula = inner_line_split[1].lstrip().rstrip()
                            fragment_dict[fragment_name]['formula'] = fragment_formula
                        #read the fragment retention time
                        if inner_read_line == outer_read_line + 2:
                            retention_time = inner_line_split[1].lstrip().rstrip()
                            retention_time = np.fromstring(retention_time,dtype=float,sep=' ')
                            fragment_dict[fragment_name]['rt'] = retention_time
                        #read the fragment list of mz values to integrate
                        if inner_read_line == outer_read_line + 3:
                            mzs_to_integrate = inner_line_split[1].lstrip().rstrip()
                            mzs_to_integrate = np.fromstring(mzs_to_integrate,dtype=float,sep=' ')
                            fragment_dict[fragment_name]['mz'] = mzs_to_integrate


    return(fragment_dict)
