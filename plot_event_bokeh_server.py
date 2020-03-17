# The call to this is:
#     bokeh serve --show plot_event_bokeh_server.py
#     the call must be made in the terminal (not in python)

import importlib
import numpy as np
import pickle #allows for saving of dictionary files
import copy
import pdb
import re
import get_directory
importlib.reload(get_directory)

from bokeh.io import curdoc
from bokeh.layouts import column, row, layout, widgetbox
from bokeh.models import ColumnDataSource, Slider, TextInput
from bokeh.plotting import figure, output_file
from bokeh.events import DoubleTap

mz_plot = ['tic','blank','blank','blank']
mz_colors = ['red','blue','green','purple']

#retrieve file directory
retrieve_directory_method = 'gui_file' #specifies you want to select the file with the gui
#    options are: 'manual', 'gui', 'manual_file', 'gui_file'
file_path = get_directory.get_directory(retrieve_directory_method)
#    returns the path to the file including the filename and extension
file_directory = re.sub('/[^/]*$','',file_path)+'/'
#    removes the filename and extension, leaving the terminating /
#        '/[^/]*$': '/' -> match '/'; ''[^/]*'' -> any character except / any number of times; '$'' -> end of string

#get the file name
regex_pattern = re.compile('/[^/]*$')
filename_regex = regex_pattern.search(file_path)
filename = filename_regex[0]
filename = re.sub('/','',filename)

sample_name = filename.split('.')[0]
input_data_file = file_directory + 'processed_data.p'
output_plot_file = file_directory + 'plot1.html'

output_file(output_plot_file, title=sample_name, mode='inline')
#    title sets the title shown on the browser tab
#    mode=inline allows it to produce the plot without being connected to the internet
#        CDN loading is the default and requires an internet connection
plot=figure(title='ion counts vs. time', x_axis_label='retention index',y_axis_label='ion counts',plot_width=950,plot_height=300)

file_object = open(input_data_file,'rb')
file_data = pickle.load(file_object)
file_object.close()

source = {}
mz_text = {}
legend = {}

source_dict_plot1 = {}
source_dict_plot2 = {}
source_dict_plot3 = {}
source_dict_plot4 = {}

x_data_source = 'ri'
#x_data_source = 'sats'
#    can be 'ri' for retention indices or 'sats' for scan acquisition times

x_data = file_data[sample_name][x_data_source] #this does not change with mz so it is set outside the loop
blank_data = np.zeros(len(x_data))
for i in range(0,len(blank_data)):
    blank_data[i] = np.nan


#change the key names to strings
source_dict = file_data[sample_name]['ics_smooth_bc'] #bc: baseline-corrected
source_dict_keys = list(source_dict.keys())
for key in source_dict_keys:
    new_key = str(key) #change the float key to a string
    new_key = re.sub('\..*$','',new_key) #remove the decimal and everything following
    source_dict[new_key] = source_dict.pop(key) #new key on left of equal sign, old key on right
source_dict['blank'] = blank_data

#set the x data and the initially displayed values
source_dict_plot1['x'] = x_data
source_dict_plot1['y'] = source_dict[mz_plot[0]]

source_dict_plot2['x'] = x_data
source_dict_plot2['y'] = source_dict[mz_plot[1]]

source_dict_plot3['x'] = x_data
source_dict_plot3['y'] = source_dict[mz_plot[2]]

source_dict_plot4['x'] = x_data
source_dict_plot4['y'] = source_dict[mz_plot[3]]


source1 = ColumnDataSource(data=source_dict_plot1) #for bokeh widgets it is stored in a ColmnDataSource object
source2 = ColumnDataSource(data=source_dict_plot2) #for bokeh widgets it is stored in a ColmnDataSource object
source3 = ColumnDataSource(data=source_dict_plot3) #for bokeh widgets it is stored in a ColmnDataSource object
source4 = ColumnDataSource(data=source_dict_plot4) #for bokeh widgets it is stored in a ColmnDataSource object

labels = ['label1','label2','label3','label4']

def update_mz_trace1(attrname, old, new):
    mz = mz_text[0].value
    y = source_dict[mz]
    source1.data = dict(x=x_data, y=y)

def update_mz_trace2(attrname, old, new):
    mz = mz_text[1].value
    y = source_dict[mz]
    source2.data = dict(x=x_data, y=y)

def update_mz_trace3(attrname, old, new):
    mz = mz_text[2].value
    y = source_dict[mz]
    source3.data = dict(x=x_data, y=y)

def update_mz_trace4(attrname, old, new):
    mz = mz_text[3].value
    y = source_dict[mz]
    source4.data = dict(x=x_data, y=y)

source_list = [source1,source2,source3,source4]
update_function_list = [update_mz_trace1,update_mz_trace2,update_mz_trace3,update_mz_trace4]

for j in [0,1,2,3]:
    plot.line('x','y',source=source_list[j],color=mz_colors[j]) #update the plot object for the current mz
    mz_text[j] = TextInput(title=mz_colors[j], value=str(mz_plot[j])) #the textbox widget, the value must be a string
    mz_text[j].on_change('value',update_function_list[j])
    #mz_text[j] = Select(title=mz_colors[j], value=str(mz_plot[j]), options=list(source_dict.keys()))

# intensity vs. mz at specified time################

#Make the ion-count vs. mz plot for each scan acquisition time
source_dict_timekeys = file_data[sample_name]['ics_smooth_timekeys']
source_dict_timekeys_keys = list(source_dict_timekeys.keys())

#set the x values to all mz values and the initial y value to the first recorded intensities for each mz
source_dict_timekeys_plot = {}
source_dict_timekeys_plot['x'] = file_data[sample_name]['mz_vals']
test_time_value = list(source_dict_timekeys.keys())[0]
source_dict_timekeys_plot['y'] = source_dict_timekeys[test_time_value]

#convert the keys into strings
source_dict_timekeys_keys = list(source_dict_timekeys.keys())
# for key in source_dict_timekeys_keys:
#     new_key = str(key) #change the float key to a string
#     #new_key = re.sub('\..*$','',new_key) #remove the decimal and everything following
#     source_dict_timekeys[new_key] = source_dict_timekeys.pop(key) #new key on left of equal sign, old key on right

source_timekeys = ColumnDataSource(data=source_dict_timekeys_plot)
plot2 = figure(title='ion counts vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=950,plot_height=300)
plot2.vbar(x='x', bottom=0, width=0.5, top='y',color='firebrick',source=source_timekeys)

#double-click callback
def callback(event):
    rt_click = event.x
    subtracting_click_time = np.array(source_dict_timekeys_keys) - rt_click
    rt_index = np.argmin(abs(subtracting_click_time))
    rt = source_dict_timekeys_keys[rt_index]
    # x_index = 'x'+'%.3f'%(rt)
    # y_index = 'y'+'%.3f'%(rt)
    x = file_data[sample_name]['mz_vals']
    y = source_dict_timekeys[rt]
    source_timekeys.data = dict(x=x, y=y)
plot.on_event(DoubleTap, callback)




# Set up layouts and add to document
text_box1 = widgetbox(mz_text[0], width=175, height=20)
text_box2 = widgetbox(mz_text[1], width=175, height=20)
text_box3 = widgetbox(mz_text[2], width=175, height=20)
text_box4 = widgetbox(mz_text[3], width=175, height=20)

l = layout([
  [text_box1,text_box2,text_box3,text_box4],
  [plot],
  [plot2],
], sizing_mode='fixed')

curdoc().add_root(l)
