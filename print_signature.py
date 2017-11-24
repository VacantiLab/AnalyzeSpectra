
import numpy as np
import importlib
from tkinter import Tk #allows for asking for a directory through a GUI
from tkinter.filedialog import askdirectory #allows for asking for a directory through a GUI
import find_closest
importlib.reload(find_closest)

#retrieve the file_directory
retrieve_directory_method = 'manual'
file_directory = get_directory.get_directory(retrieve_directory_method)

#initialize the filename and retention index for testing the function
filename = 'tbdms01_t47d_wt.CDF'
sample_name = 'tbdms01_t47d_wt'
input_data_file = file_directory + 'processed_data.p'
ri = 1393
ri_array = file_data[sample_name]['ri']

#open the data specified by the filename
file_object = open(input_data_file,'rb')
file_data = pickle.load(file_object)
file_object.close()

#find the closest recorded retention index to that specified
ri_closest = find_closest.find_closest(ri,ri_array)[1]

#initialize coelution arrays for testing the function
coelution_array = file_data[sample_name]['coelution_dictionary'][ri_closest]
coelution_val = file_data[sample_name]['coelution_dicionary_values'][ri_closest]

#initialize the array containing the peaks you considered when looking for groups
peaks_considered = np.array([])

#initialize the vector containing the mz values that belong in groups and the corresponding ion counts at those values
#    groups dictionary: each group will be named by the initial value and members will be all mz values that are part of that group
#    values dictionary: the key names correspond to the groups dictionary, values are corresponding ion counts
groups = dict()
values = dict()

#iterate through the coleution array
j = 0
for mz in coelution_array:
    #determined if you already considered whether the current mz is a member of a group
    peak_considered = mz in peaks_considered
    #if you have not already considered it, consider it now
    if not peak_considered:
        peaks_considered = np.append(peaks_considered,mz) #add the mz value to the list of considered values
        peak_count = 0 #initialize the peak count within this group
        peak_present = True #a peak is present for this mz value because it came from the coelution array
        #keep checking to see if the next mz value is in the coeluting peaks until you find it is not
        while peak_present == True:
            peak_count = peak_count+1 #increase the peak count, because you entered this loop
            peaks_considered = np.append(peaks_considered,mz+peak_count) #you just considered another peak, include it in the list
            #if you have three consecutive mz values eluting, you have a group
            if peak_count==3:
                groups[mz] = np.array([mz,mz+1,mz+2]) #set the initial group array stored in the group dictionary
                values[mz] = np.array([coelution_val[j],coelution_val[j+1],coelution_val[j+2]]) #store the corresponding ion count values for each mz in the group
            #continue adding to that group if it grows
            if peak_count > 3:
                groups[mz] = np.append(groups[mz],mz+peak_count-1)
                values[mz] = np.append(values[mz],coelution_val[j+peak_count-1])
            #determine if the next peak is present and if you should go through the while loop again
            peak_present = mz+peak_count in coelution_array
    j = j+1

#determine the mz values that start a group of mz values
mz_groups = list(dict.keys(groups))

#calculate the signature array
#    this array will have the initial mz value of a group followed by the total ion count of that group, followed by the relative abundance of each member
#    for example a signature array of 175 203456 0.63 0.25 0.1 0.02 233 1407893 0.5 0.4 0.1 means there are 2 groups: 175 176 177 178 and 233 234 235
#        for 175: there are 203456 total ion counts and the corresponding relative abundances of M0 to M3 are 0.63 0.25 0.1 0.02
signature_array = np.array([])
for mz in mz_groups:
    signature_array = np.append(signature_array,mz)
    values_array = values[mz]
    group_tic = np.sum(values_array) #get the total sum intensity of all group members
    signature_array = np.append(signature_array,group_tic)
    values_array_norm = values_array/np.sum(values_array)
    indices_to_iterate = np.arange(0,len(values_array),1) #start at 1 because you do not want the 0 index here
    for i in indices_to_iterate:
        signature_array = np.append(signature_array,values_array_norm[i])

#convert items in the signature_array to strings for printing into a .txt file
sig_array_str = np.array([])
for item in signature_array:
    #round values less than one to 3 decimal places
    if item < 1:
        item = np.round(item,3)
    sig_array_str = np.append(sig_array_str,item.astype('str'))

#write to the output file
file_path = file_directory + 'fragment_signature.txt'
with open(file_path,'w') as sig_file:
    j = 0
    for item in sig_array_str:
        j = j+1
        if j < len(sig_array_str):
            sig_file.write(item)
            sig_file.write(' ')
        if j == len(sig_array_str):
            sig_file.write(item)
