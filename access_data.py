import importlib #allows fresh importing of modules
import pdb #python debugger
#import plotly #for plotting, not used here?
from os import listdir #needed to get a list of files
import re #for regular expressions
import bokeh.plotting as bkp #allows for making interactive plots with bokeh
import numpy as np #this is numpy, allows for data frame and matrix handling
import pandas #a module which allows for making data frames
import copy
import pickle #allows for saving of dictionary files
from tkinter import Tk #allows for asking for a directory through a GUI
from tkinter.filedialog import askdirectory #allows for asking for a directory through a GUI
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
import find_ri
importlib.reload(find_ri)
import get_ri_keys_dict
importlib.reload(get_ri_keys_dict)

#ask for the directory where the netCDF and library.txt files are
root = Tk()
root.withdraw() #closes the tkinter GUI window because the rest of the program is not run through the GUI
file_directory = askdirectory() + '/'
root.update() #required so the directory request dialog box disappears and does not freeze

#Get a list of all of the files in the specified directory
files = listdir(file_directory)

library_processed = 'library.p' in files

if not library_processed:
    #process the library
    print('processing library...')
    fragment_dict,fragment_list = fragment_library.fragment_library()

    #Save the library into a python readable file
    output_library_file = file_directory + 'library.p'
    with open(output_library_file,'wb') as library_file_object:
        pickle.dump(fragment_dict,library_file_object)

if library_processed:
    #load the processed library
    print('opening previously processed library...')
    input_library_file = file_directory + 'library.p'
    with open(input_library_file,'rb') as library_file_object:
        fragment_dict = pickle.load(library_file_object)
        fragment_list = list(dict.keys(fragment_dict))

#Remove filenames that are not NetCDF files
netcdf_pattern = re.compile('.cdf$|.netcdf$',re.IGNORECASE)
    #creates a regex pattern matching '.cdf' or '.netcdf' at the end of a string (specified by $), where the case is not important

filename_holder = copy.copy(files)
for filename in filename_holder:
    is_netcdf = bool(re.search(netcdf_pattern,filename)) #determines if the pattern is in the filename
    if not is_netcdf:
        files.remove(filename)

#Create the folder for outputting results
output_plot_directory,output_directory = create_output_directory.create_output_directory()

#Initialize a dictionary to contain all of the outputs of integrating each NetCDF file in the specified directory
file_data = {}

#Iterate through each NetCDF file and process the data
for filename in files:
    print(filename + ':')
    file_path = file_directory+filename

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
    peak_start_i_dict,peak_end_i_dict,x_data_numpy,p) = process_ms_data.process_ms_data(sat,ic_df,output_plot_directory,n_scns,mz_vals)
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

    #transform scan acquisition times to retention indices
    if filename == 'alkanes.CDF':
        ri_array = find_ri.find_ri(ic_smooth_dict,mz_vals,sat)

    #inverty the ic_smooth_dict so that retention indices are the keys and a vector of intensities for each mz are the items
    ic_smooth_dict_timekeys = get_ri_keys_dict.get_ri_keys_dict(ic_smooth_dict,ri_array,mz_vals)

    #integrate fragments in library
    print('    integrating fragment mass isotopomers listed in library...')
    fragment_dict_complete = integrate_peaks.integrate_peaks(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,
                                                    peak_start_i_dict,peak_end_i_dict,x_data_numpy,fragment_dict,fragment_list)
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
    file_data[filename] = {}
    file_data[filename]['fragments'] = copy.deepcopy(fragment_dict_complete)
    file_data[filename]['ics_smooth_bc'] = ic_smooth_dict #bc stands for baseline-corrected
    file_data[filename]['ics_smooth_timekeys'] = ic_smooth_dict_timekeys
    file_data[filename]['sats'] = x_data_numpy
    file_data[filename]['ri'] = ri_array
    file_data[filename]['mz_vals'] = mz_vals
    file_data[filename]['peak_beginnings'] = peak_start_t_dict
    file_data[filename]['peak_endings'] = peak_end_t_dict

#Save the output data into a python readable file
output_data_file = file_directory + 'processed_data.p'
file_object = open(output_data_file,'wb')
pickle.dump(file_data,file_object)
file_object.close()

#Print output to a text file_data
print_integrated_peaks.print_integrated_peaks(file_directory,files,fragment_list,file_data)

print('Data processed successfully.')
