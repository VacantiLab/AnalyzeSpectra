
import numpy as np
from tkinter import Tk #allows for asking for a directory through a GUI
from tkinter.filedialog import askdirectory #allows for asking for a directory through a GUI

coelution_array = np.array([174,175,176,177,178,180,192,233,234,235,236,303,307])
coelution_val = np.array([1003,506,102,23,500,607,3067,1001,200,43,76,89,87])

j=0
#initialize the array containing the peaks you searched for within the coelution array
peaks_considered = np.array([])

#initialize the vector containing the mz values that belong in groups
groups = dict()
values = dict()

#iterate through the coleution array
j = 0
for mz in coelution_array:
    #determined if you already searched for the current mz
    peak_considered = mz in peaks_considered
    #if you have not already considered it, consider it now
    if not peak_considered:
        peaks_considered = np.append(peaks_considered,mz) #add the mz value to the list of considered values
        peak_count = 0 #initialize the peak count within this group
        peak_present = True #a peak is present for this mz value because it came from the coelution array
        while peak_present == True:
            peak_count = peak_count+1 #increase the peak count, because you entered this loop
            peaks_considered = np.append(peaks_considered,mz+peak_count) #you just considered another peak, include it in the list
            #if you have three peaks in a row, you have a group
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

mz_groups = list(dict.keys(groups))
signature_array = np.array([])

for mz in mz_groups:
    signature_array = np.append(signature_array,mz)
    values_array = values[mz]
    values_array_norm = values_array/np.sum(values_array)
    for i in np.arange(1,len(values_array),1):
        signature_array = np.append(signature_array,values_array_norm[i])

#ask for the directory where to print the fragment signature
root = Tk()
root.withdraw() #closes the tkinter GUI window because the rest of the program is not run through the GUI
file_directory = askdirectory() + '/'
root.update() #required so the directory request dialog box disappears and does not freeze

file_path = file_directory + 'fragment_signature.txt'
with open(file_path,'w') as sig_file:
    for item in signature_array:
        sig_file.write(item)
        sig_file.write('\t')
