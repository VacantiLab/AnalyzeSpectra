def integrate(corrected=True):

    #This script expects two tab-delimited text files to be stored where the data is stored
    #    The first is files_to_batch.txt
    #        This file has two columns, the first "file_name" and the second "batch"
    #        The file_name is the name of the file including the extension
    #        The batch is a number indicating which files will be analyzed with the same alkane file
    #    The second is "batch_to_alkane.txt"
    #        This file has two columns, the first "batch" and the second "alkane_file"
    #        The "batch" is the number of the batch
    #        The "alkane_file" is the alkane file corresponding to that batch number
    #This script can be run without either text file described above
    #    In that case it would reference only a sinlge alkane file that must be named: alkanes.CDF

    import importlib #allows fresh importing of modules
    from pdb import set_trace #python debugger
    #import plotly #for plotting, not used here?
    from os import listdir #needed to get a list of files
    import re #for regular expressions
    import bokeh.plotting as bkp #allows for making interactive plots with bokeh
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    import copy
    import pickle #allows for saving of dictionary files
    import organize_ms_data
    importlib.reload(organize_ms_data)
    import process_ms_data
    importlib.reload(process_ms_data)
    import create_output_directory
    importlib.reload(create_output_directory)
    import integrate_peaks
    importlib.reload(integrate_peaks)
    import BinData
    importlib.reload(BinData)
    import fragment_library
    importlib.reload(fragment_library)
    import print_integrated_peaks
    importlib.reload(print_integrated_peaks)
    import find_ri_conversion
    importlib.reload(find_ri_conversion)
    import get_ri_keys_dict
    importlib.reload(get_ri_keys_dict)
    import calc_coelut
    importlib.reload(calc_coelut)
    import convert_rt_ri
    importlib.reload(convert_rt_ri)
    import get_directory
    importlib.reload(get_directory)
    import locate_overlap
    importlib.reload(locate_overlap)
    import GetFileBatch
    importlib.reload(GetFileBatch)

    #retrieve file directory
    retrieve_directory_method = 'gui'
    file_directory = get_directory.get_directory(retrieve_directory_method)

    #Get a list of all of the files in the specified directory
    files = listdir(file_directory)

    library_processed = 'library.p' in files

    if not library_processed:
        #process the library
        print('processing library...')
        metabolite_dict,metabolite_list = fragment_library.fragment_library(file_directory)

    if library_processed:
        #load the processed library
        print('opening previously processed library...')
        input_library_file = file_directory + 'library.p'
        with open(input_library_file,'rb') as library_file_object:
            metabolite_dict = pickle.load(library_file_object)
            metabolite_list = list(dict.keys(metabolite_dict))
        #Check for new metabolites
        metabolite_dict,metabolite_list = fragment_library.fragment_library(file_directory,metabolite_dict)

    #Save the library into a python readable file
    output_library_file = file_directory + 'library.p'
    with open(output_library_file,'wb') as library_file_object:
        pickle.dump(metabolite_dict,library_file_object)


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

    #get the mappings of the sample file names to the propper batch and the batch names to the propper alkane file
    file_name_to_batch,batch_to_alkane,alkane_files,batches = GetFileBatch.GetFileBatch(file_directory)
    #    This file expects two tab-delimited text files to be stored where the data is stored
    #        The first is files_to_batch.txt
    #            This file has two columns, the first "file_name" and the second "batch"
    #            The file_name is the name of the file including the extension
    #            The batch is a number indicating which files will be analyzed with the same alkane file
    #        The second is "batch_to_alkane.txt"
    #            This file has two columns, the first "batch" and the second "alkane_file"
    #            The "batch" is the number of the batch
    #            The "alkane_file" is the alkane file corresponding to that batch number

    #remove the alkane files from the file list
    for alkane_file in alkane_files:
        files.remove(alkane_file)

    #Create the folder for outputting results
    output_plot_directory,output_directory = create_output_directory.create_output_directory()

    #Initialize a dictionary to contain all of the outputs of integrating each NetCDF file in the specified directory
    file_data = {}

    #Initialize an array to contain all sample names
    samples_all = np.array([])

    #Separate the file names by their batches and integrate them using the propper alkane file
    #iterate through the batches and analyze
    for batch in batches:

        #Get a list of the filenames in the current batch
        indices = file_name_to_batch['batch']==batch
        files = file_name_to_batch['file_name'][indices]
        files = np.array(files)

        #Get the name of the corresponding alkane file
        #    the variable batches becomes a different data type if there is only one batch
        #        thus there are the if statements to handle it correctly
        if len(batches)>1:
            alkane_file_index = batch_to_alkane['batch']==batch
        if len(batches)==1:
            alkane_file_index = 0
        alkane_file = batch_to_alkane['alkane_file'][alkane_file_index]
        if len(batches)>1:
            alkane_file = np.array(alkane_file) #it is converted to a np array so the next line works
            alkane_name = alkane_file[0].split('.')[0] #removes the .CDF from the end of the filename
        alkane_name = alkane_file.split('.')[0] #removes the .CDF from the end of the filename

        #Add the name of the alkane file to the beginning of the list
        files = np.insert(files,0,alkane_file)

        #Iterate through each NetCDF file and process the data
        i=0
        samples = copy.copy(files) #the sample names will be the filenames without the extensions
        #be sure to start with the alkane file
        for filename in files:
            print(filename + ':')
            file_path = file_directory+filename

            #get the sample names and store them
            sample_name = filename.split('.')[0]
            samples[i] = sample_name

            print('    accessing and organizing m/z, scan acquisition time, and ion count data...')
            ic_df,sat,n_scns,mz_vals,tic = organize_ms_data.organize_ms_data(file_path)
                #ic_df: ion count data frame
                #sat: scan acquisition times
                #n_scns: number of scans
                #mz_vals: the m/z values scanned

            #bin the data as specified: move to organize_ms_data
            ic_df = BinData.BinData(ic_df,1.0005)
                #the second input is the width of the bin

            #reassign the mz values due to the binning
            mz_vals = np.sort(np.array(list(ic_df.index.values)))

            #Process ms data
            print('    subtracting baselines and smoothing...')
            (ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,
            peak_start_i_dict,peak_end_i_dict,x_data_numpy,peak_i_dict,
            peak_max_dict,p,peak_sat_dict,ic_dict) = process_ms_data.process_ms_data(sat,ic_df,output_plot_directory,n_scns,mz_vals)
                 #ic_smooth_dict: a dictionary containing the smoothed and baseline corrected ion count data for each m/z value
                 #peak_start_t_dict: a dictionary with all of the peak beginning times for each m/z ion count plot
                 #peak_end_t_dict: a dictionary with all of the peak ending times for each m/z ion count plot
                 #peak_start_i_dict: a dictionary with all of the peak beginning indices for each m/z ion count plot
                 #peak_end_i_dict: a dictionary with all of the peak ending indices for each m/z ion count plot
                 #x_data_numpy: scan acquisition time values
                 #p: the plot (bokeh) object

            #add the total ion count to the dictionary
            #    note it is not smoothed because it does not really make sense to smooth total ion count data
            ic_smooth_dict['tic'] = tic
            ic_dict['tic'] = tic

            #calculate peak overlap dictionary
            print('    finding coeluting peaks ...')
            peak_overlap_dictionary = locate_overlap.locate_overlap(ic_smooth_dict,peak_start_i_dict,peak_end_i_dict,mz_vals,peak_max_dict)

            #the first sample must always be alkanes - plan to make this optional later
            #find the retention time to retention index conversion
            #    one array is retention indices and the other is corresponding retention times
            if sample_name == alkane_name:
                #calculate the coelution dictionary with the scan acquisition times as keys
                #    coelution_dict_sat has keys of sat's and arrays of mz's whoe peaks elute at those sat's
                #    coelution_dict_val is the same except the arrays are the corresponding intensity values of the eluting peaks at the sat of the key
                coelut_dict_sat,coelut_dict_val_sat = calc_coelut.calc_coelut(peak_sat_dict,mz_vals,sat,ic_smooth_dict,peak_overlap_dictionary)
                ri_sat,ri_rec = find_ri_conversion.find_ri_conversion(ic_smooth_dict,mz_vals,sat,coelut_dict_sat,coelut_dict_val_sat,sample_name)

            #convert the retention times of the current sample to retention indices
            #    doing this for each sample allows for samples with differing quantities of scan acquisition times
            #    to be analyzed with the same alkane sample for retention index calculation
            ri_array = convert_rt_ri.convert_rt_ri(ri_sat,ri_rec,sat)

            #find ri of each peak
            peak_ri_dict = dict()
            peak_start_ri_dict = dict()
            peak_end_ri_dict = dict()
            for mz_val in mz_vals:
                peak_loc_ind = peak_i_dict[mz_val]
                peak_start_ind = peak_start_i_dict[mz_val]
                peak_end_ind = peak_end_i_dict[mz_val]
                peak_ri_dict[mz_val] = ri_array[peak_loc_ind]
                peak_start_ri_dict[mz_val] = ri_array[peak_start_ind]
                peak_end_ri_dict[mz_val] = ri_array[peak_end_ind]

            #invert the ic_smooth_dict so that retention indices are the keys and a vector of intensities for each mz are the items
            ic_smooth_dict_timekeys = get_ri_keys_dict.get_ri_keys_dict(ic_smooth_dict,ri_array,mz_vals)

            #calculate the coelution dictionary with the retention indices as keys
            #    coelution_dict has keys of ri's and arrays of mz's whoe peaks elute at those ri's
            #    coelution_dict_val is the same except the arrays are the corresponding intensity values of the eluting peaks at the ri of the key
            coelut_dict,coelut_dict_val = calc_coelut.calc_coelut(peak_ri_dict,mz_vals,ri_array,ic_smooth_dict,peak_overlap_dictionary)

            #integrate fragments in library
            print('    integrating fragment mass isotopomers listed in library...')
            metabolite_dict_complete = integrate_peaks.integrate_peaks(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,
                                                            peak_start_i_dict,peak_end_i_dict,x_data_numpy,metabolite_dict,
                                                            metabolite_list,ri_array,mz_vals,coelut_dict,coelut_dict_val,sample_name)
            #fragment_dict: a dictionary containing information (including the mass isotopomer distributions) of each integrated metabolite fragment

            #Store the processed data for each filename in a dictionary
            #This is a dictionary of dicionaries tree with the following structure:
            #    file_data
            #    ....filename (there is one for each file integrated)
            #        ....fragments
            #            ....fragment_id (there is one for every fragment integrated)
            #                ....formula
            #                ....rt
            #                ....mz
            #                ....areas
            #                    ....mz(there is one for every mz integrated for the fragment_id)
            #                ....mid
            #                    ....mz(there is one for every mz integrated for the fragment_id)
            #                ....total_area
            #        ....ics_smooth_bc
            #            ....mz_value (there is one of these for every mz scanned)
            #        ....sats
            #        ....peak_beginnings
            #            ....mz_value (there is one of these for every mz scanned)
            #        ....peak_endings
            #            ....mz_value (there is one of these for every mz scanned)
            file_data[sample_name] = {}
            file_data[sample_name]['metabolites'] = copy.deepcopy(metabolite_dict_complete)
            file_data[sample_name]['ics_smooth_bc'] = ic_smooth_dict #bc stands for baseline-corrected
            file_data[sample_name]['ics_smooth_timekeys'] = ic_smooth_dict_timekeys
            file_data[sample_name]['ics'] = ic_dict
            file_data[sample_name]['sats'] = x_data_numpy
            file_data[sample_name]['ri'] = ri_array
            file_data[sample_name]['mz_vals'] = mz_vals
            file_data[sample_name]['peak_beginnings'] = peak_start_t_dict
            file_data[sample_name]['peak_endings'] = peak_end_t_dict
            file_data[sample_name]['peak_indices'] = peak_i_dict
            file_data[sample_name]['peak_ris'] = peak_ri_dict
            file_data[sample_name]['peak_maxes'] = peak_max_dict
            file_data[sample_name]['peak_start_ris'] = peak_start_ri_dict
            file_data[sample_name]['peak_end_ris'] = peak_end_ri_dict
            file_data[sample_name]['coelution_dictionary'] = coelut_dict
            file_data[sample_name]['coelution_dicionary_values'] = coelut_dict_val
            file_data[sample_name]['peak_overlap_dictionary'] = peak_overlap_dictionary

            i=i+1

        #Remove the alkane file from the sample list and update the samples_all list
        samples = np.delete(samples,0)
        samples_all = np.append(samples_all,samples)

    #Save the output data into a python readable file
    output_data_file = file_directory + 'processed_data.p'
    file_object = open(output_data_file,'wb')
    pickle.dump(file_data,file_object)
    file_object.close()

    #Print output to a text file_data
    print_integrated_peaks.print_integrated_peaks(file_directory,samples_all,metabolite_list,file_data,corrected)

    print('Data processed successfully.')

    return()
