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
    #    metabolite: pyruvate
    #    retention index: 1223
    #    peak profile: 151 152 174 175
    #    fragment: pyr174
    #    formula: C6H12N1O3Si1
    #    mzs to integrate: 174 175 176 177 178
    #    ...
    with open(file_name_read, 'r') as read_file:
        #read through the lines in the file one by one
        file_line_n = 0
        for line in read_file:
            file_line_n = file_line_n + 1 #iterate the line number
            line_split = line.split(':') #split the line into a list of strings, colons mark separations
            #if the first word on the line is 'metabolite', gather the fragment information
            if line_split[0]=='metabolite':
                before_blankline_border = True
                metabolite_name = line_split[1].lstrip().rstrip() #remove white space characters from the left and right of the fragment name
                print('    ' + metabolite_name)
                metabolite_dict[metabolite_name] = dict() #initialize a dictionary for the current metabolite
                metabolite_dict[metabolite_name]['fragments'] = dict() #initialize a dictionary for the fragments of the current metabolite
                metabolite_dict[metabolite_name]['peak_profile'] = np.array([]) #initialize an array for the peak profile of the current metabolite
                #read through the lines in the file one by one again, stopping once you are passed where you stopped previously
                #    this puts you one past the metabolite line
                with open(file_name_read, 'r') as metabolite_read_file:
                    metabolite_line_n = 0
                    for metabolite_line in metabolite_read_file:
                        metabolite_line_n = metabolite_line_n + 1
                        metabolite_line_split = metabolite_line.split(':')
                        metabolite_line_title = metabolite_line_split[0].lstrip().rstrip()
                        #if you are one passed the metabolite line, you are at the retention index line
                        if metabolite_line_n == file_line_n + 1:
                            metabolite_line_item = metabolite_line_split[1].lstrip().rstrip() #must be defined within the if statement because it causes an error for an empty line
                            retention_index = metabolite_line_item
                            retention_index = np.fromstring(retention_index,dtype=float,sep=' ')
                            metabolite_dict[metabolite_name]['ri'] = retention_index
                        #if you are two passed the metabolite line, you are at the peak profile line
                        if metabolite_line_n == file_line_n + 2:
                            metabolite_line_item = metabolite_line_split[1].lstrip().rstrip()
                            peak_profile = metabolite_line_item
                            peak_profile = np.fromstring(peak_profile,dtype=float,sep=' ')
                            metabolite_dict[metabolite_name]['peak_profile'] = peak_profile

                        #if you run into an empty line after the line marking the current metabolite, take note
                        #    you should not proceed with collecting fragment information because you are on the next metabolite
                        if (metabolite_line in ['\n', '\r\n']) & (metabolite_line_n > file_line_n):
                            before_blankline_border = False

                        #if the line begins with 'fragment', this is the name of a fragment of the metabolite that is integrated for MID calculation
                        #    if this is the case, read through the file again to record the fragment information
                        #    you must not have passed a blank space line prior to seeing the fragment indicator, elsewise you are on a new metabolite
                        #    you must be passed the current metabolite indicator line, elsewise you are on the previous metabolite
                        if (metabolite_line_title == 'fragment') & before_blankline_border & (metabolite_line_n > file_line_n):
                            with open(file_name_read, 'r') as fragment_read_file:
                                fragment_line_n = 0
                                for fragment_line in fragment_read_file:
                                    fragment_line_n = fragment_line_n + 1
                                    fragment_line_split = fragment_line.split(':')
                                    #once you have the fragment name, you can initialize a dictionary containing all of the fragment information
                                    #    you can also initialize the numpy arrays that are members of that dictionary
                                    if fragment_line_n == metabolite_line_n:
                                        fragment_line_item = fragment_line_split[1].lstrip().rstrip() #must be defined within the if statement because it causes an error for an empty line
                                        fragment_name = fragment_line_item
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name] = dict()
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['areas'] = np.array([])
                                    #when you are at the formula line, record the formula and calculate the inverted correction matrix
                                    if fragment_line_n == metabolite_line_n + 1:
                                        fragment_line_item = fragment_line_split[1].lstrip().rstrip()
                                        frag_formula = fragment_line_item
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['formula'] = frag_formula
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['CM_i'] = create_correction_matrix.create_correction_matrix(frag_formula)
                                    #when at the appropriate line, record the mz's that will be integrated
                                    if fragment_line_n == metabolite_line_n + 2:
                                        fragment_line_item = fragment_line_split[1].lstrip().rstrip()
                                        mzs_to_integrate = np.fromstring(fragment_line_item,dtype=float,sep=' ')
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['mzs_to_integrate'] = mzs_to_integrate
                                    if fragment_line_n == metabolite_line_n + 3:
                                        fragment_line_item = fragment_line_split[1].lstrip().rstrip()
                                        metabolite_atoms = fragment_line_item
                                        metabolite_dict[metabolite_name]['fragments'][fragment_name]['metabolite_atoms'] = metabolite_atoms


    #get a list of all of the metabolites
    metabolite_list = list(dict.keys(metabolite_dict))

    for metabolite_name in metabolite_list:
        #calculate the natural mass isotopomer distrubutions for each fragment
        fragment_list = list(dict.keys(metabolite_dict[metabolite_name]['fragments']))
        for z in fragment_list:
            metabolite_dict[metabolite_name]['fragments'][z]['natural_mid'] = calc_natural_mid.calc_natural_mid(metabolite_dict[metabolite_name]['fragments'][z]['formula'])

    return(metabolite_dict,metabolite_list)
