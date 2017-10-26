def fragment_library():

    import importlib #allows fresh importing of modules
    import pdb
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas
    import re
    import calc_natural_mid
    importlib.reload(calc_natural_mid)
    import create_correction_matrix
    importlib.reload(create_correction_matrix)

    #file_name_read = '/Users/Nate/Desktop/netcdf_test/tbdms_lib.txt'
    file_name_read = '/Users/Nate/Dropbox/Research/Lehtio_Laboratory/Projects/metabolite_integration_tool/netcdf_test/library.txt'
    metabolite_dict = dict()

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
        file_line_n = 0
        for line in read_file:
            file_line_n = file_line_n + 1 #iterate the line number
            line_split = line.split(':') #split the line into a list of strings, colons mark separations
            #if the first word on the line is 'fragment', gather the fragment information
            if line_split[0]=='metabolite':
                metabolite_name = line_split[1].lstrip().rstrip() #remove white space characters from the left and right of the fragment name
                print('    ' + metabolite_name)
                metabolite_dict[metabolite_name] = dict() #initialize a dictionary for the current fragment
                metabolite_dict[metabolite_name]['fragments'] = dict()
                metabolite_dict[metabolite_name]['peak_profile'] = np.array([])
                with open(file_name_read, 'r') as metabolite_read_file:
                    metabolite_line_n = 0
                    for metabolite_line in metabolite_read_file:
                        metabolite_line_n = metabolite_line_n + 1
                        metabolite_line_split = metabolite_line.split(':')
                        metabolite_line_title = metabolite_line_split[0].lstrip().rstrip()
                        metabolite_line_item = metabolite_line_split[1].lstrip().rstrip()
                        if metabolite_line_n == file_line_n + 1:
                            retention_index = metabolite_line_item
                            retention_index = np.fromstring(retention_index,dtype=float,sep=' ')
                            metabolite_dict[metabolite_name]['ri'] = retention_index
                        if metabolite_line_n == file_line_n + 2:
                            peak_profile = metabolite_line_item
                            peak_profile = np.fromstring(peak_profile,dtype=float,sep=' ')
                            metabolite_dict[metabolite_name]['peak_profile'] = peak_profile
                        if metabolite_line_title == 'fragment':
                            with open(file_name_read, 'r') as fragment_read_file:
                                fragment_line_n = 0
                                for fragment_line in fragment_read_file:
                                    fragment_line_n = fragment_line_n + 1
                                    fragment_line_split = fragment_line.split(':')
                                    fragment_line_item = fragment_line_split[1].lstrip().rstrip()
                                    if fragment_line_n == metabolite_line_n:
                                        fragment_name = fragment_line_item
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name] = dict()
                                    if fragment_line_n == metabolite_line_n + 1:
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['formula'] = fragment_line_item
                                    if fragment_line_n == metabolite_line_n + 2:
                                        mzs_to_integrate = np.fromstring(fragment_line_item,dtype=float,sep=' ')
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['mzs_to_integrate'] = mzs_to_integrate

    metabolite_list = list(metabolite_dict.keys())
    
    #calculate the natural mass isotopomer distrubutions for each fragment
    #fragment_list = list(metabolite_dict[metabolite_name]['fragments'].keys())
    #for z in fragment_list:
    #    metabolite_dict[metabolite_name][z]['natural_mid'] = calc_natural_mid.calc_natural_mid(metabolite_dict[metabolite_name]['fragments'][z]['formula'])

    fragment_list = np.array([])
    return(metabolite_dict,fragment_list)
