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
import organize_ms_data
importlib.reload(organize_ms_data)
import process_ms_data
importlib.reload(process_ms_data)
import create_output_directory
importlib.reload(create_output_directory)
import integrate_peaks
importlib.reload(integrate_peaks)

#Get a list of all of the NetCDF files in the specified directory
file_directory = '/Users/nate/Desktop/netcdf_test/'
files = listdir(file_directory)

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
    file_path = file_directory+filename

    ic_df,sat,n_scns,mz_vals = organize_ms_data.organize_ms_data(file_path)
        #ic_df: ion count data frame
        #sat: scan acquisition times
        #n_scns: number of scans
        #mz_vals: the m/z values scanned

    #Process ms data
    (ic_smooth_dict,ic_smooth_dict_timekeys,peak_start_t_dict,peak_end_t_dict,
    peak_start_i_dict,peak_end_i_dict,x_data_numpy,p) = process_ms_data.process_ms_data(sat,ic_df,output_plot_directory,n_scns,mz_vals)
         #ic_smooth_dict: a dictionary containing the smoothed and baseline corrected ion count data for each m/z value
         #peak_start_t_dict: a dictionary with all of the peak beginning times for each m/z ion count plot
         #peak_end_t_dict: a dictionary with all of the peak ending times for each m/z ion count plot
         #peak_start_i_dict: a dictionary with all of the peak beginning indices for each m/z ion count plot
         #peak_end_i_dict: a dictionary with all of the peak ending indices for each m/z ion count plot
         #x_data_numpy: scan acquisition time values
         #p: the plot (bokeh) object

    #integrate fragments in library
    fragment_dict = integrate_peaks.integrate_peaks(ic_smooth_dict,peak_start_t_dict,peak_end_t_dict,
                                                peak_start_i_dict,peak_end_i_dict,x_data_numpy)
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
    file_data[filename]['fragments'] = fragment_dict
    file_data[filename]['ics_smooth_bc'] = ic_smooth_dict #bc stands for baseline-corrected
    file_data[filename]['ics_smooth_timekeys'] = ic_smooth_dict_timekeys
    file_data[filename]['sats'] = x_data_numpy
    file_data[filename]['mz_vals'] = mz_vals
    file_data[filename]['peak_beginnings'] = peak_start_t_dict
    file_data[filename]['peak_endings'] = peak_end_t_dict

#Save the output data into a python readable file
output_data_file = file_directory + 'processed_data.p'
file_object = open(output_data_file,'wb')
pickle.dump(file_data,file_object)
file_object.close()

file_object = open(output_data_file,'rb')
file_data2 = pickle.load(file_object)
file_object.close()

    #show the plot
    #bkp.show(p)
