import importlib
import numpy as np
import pickle #allows for saving of dictionary files
import copy
import pdb
import re
import get_directory
importlib.reload(get_directory)

from bokeh.layouts import column, row, widgetbox, layout
from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput, Select
from bokeh.plotting import Figure, output_file, show, reset_output
from bokeh.io import curdoc

mz_plot = ['tic','blank','blank','blank']
mz_colors = ['red','blue','green','purple']

#retrieve file directory
retrieve_directory_method = 'manual'
file_directory = get_directory.get_directory(retrieve_directory_method)

filename = 'P01_SUM149_NT_siRNA.CDF'
sample_name = filename.split('.')[0]
input_data_file = file_directory + 'processed_data.p'
output_plot_file = file_directory + 'plot1.html'

output_file(output_plot_file, title=sample_name, mode='inline')
#    title sets the title shown on the browser tab
#    mode=inline allows it to produce the plot without being connected to the internet
#        CDN loading is the default and requires an internet connection
plot=Figure(title='ion counts vs. time', x_axis_label='retention index',y_axis_label='ion counts',plot_width=950,plot_height=300)

file_object = open(input_data_file,'rb')
file_data = pickle.load(file_object)
file_object.close()

source = {}
mz_text = {}
legend = {}

x_data = file_data[sample_name]['ri'] #this does not change with mz so it is set outside the loop
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

#set the x data and the initially displayed values
source_dict['x'] = x_data
source_dict['blank'] = blank_data
source_dict['y0'] = source_dict[mz_plot[0]] #the ics for the current mz
source_dict['y1'] = source_dict[mz_plot[1]]
source_dict['y2'] = source_dict[mz_plot[2]]
source_dict['y3'] = source_dict[mz_plot[3]]

#set the legend labels of the initially displayed values
#    there must be a label specified for each point, thus the mz lable is repeated for as many points as there are
source_dict['label0'] = np.repeat(mz_plot[0],len(source_dict[mz_plot[0]]))
source_dict['label1'] = np.repeat(mz_plot[1],len(source_dict[mz_plot[0]]))
source_dict['label2'] = np.repeat(mz_plot[2],len(source_dict[mz_plot[0]]))
source_dict['label3'] = np.repeat(mz_plot[3],len(source_dict[mz_plot[0]]))

source = ColumnDataSource(data=source_dict) #for bokeh widgets it is stored in a ColmnDataSource object

y = ['y0','y1','y2','y3']
labels = ['label0','label1','label2','label3']

for j in [0,1,2,3]:
    plot.line('x',y[j],source=source,color=mz_colors[j],legend=labels[j]) #update the plot object for the current mz
    mz_text[j] = TextInput(title=mz_colors[j], value=str(mz_plot[j])) #the textbox widget, the value must be a string
    #mz_text[j] = Select(title=mz_colors[j], value=str(mz_plot[j]), options=list(source_dict.keys()))

callback0 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y0 = data['y0']
    label0 = data['label0']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y0[i] = to_change_to[i]
        label0[i] = f
    }
    data['label0'] = label0
    source.change.emit();
""")

callback1 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y1 = data['y1']
    label1 = data['label1']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y1[i] = to_change_to[i]
        label1[i] = f
    }
    data['label1'] = label1
    source.change.emit();
""")

callback2 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y2 = data['y2']
    label2 = data['label2']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y2[i] = to_change_to[i]
        label2[i] = f
    }
    data['label2'] = label2
    source.change.emit();
""")

callback3 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y3 = data['y3']
    label3 = data['label3']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y3[i] = to_change_to[i]
        label3[i] = f
    }
    data['label3'] = label3
    source.change.emit();
""")

#inform bokeh to update the plot based on textbox values
mz_text[0].js_on_change('value', callback0)
mz_text[1].js_on_change('value', callback1)
mz_text[2].js_on_change('value', callback2)
mz_text[3].js_on_change('value', callback3)

#Make the ion-count vs. mz plot for each scan acquisition time
source_dict_timekeys = file_data[sample_name]['ics_smooth_timekeys']
source_dict_timekeys_keys = list(source_dict_timekeys.keys())

#set the x values to all mz values and the initial y value to the first recorded intensities for each mz
source_dict_timekeys['x'] = file_data[sample_name]['mz_vals']
test_time_value = list(source_dict_timekeys.keys())[0]
source_dict_timekeys['y'] = source_dict_timekeys[test_time_value]

#convert the keys into strings
source_dict_timekeys_keys = list(source_dict_timekeys.keys())
for key in source_dict_timekeys_keys:
    new_key = str(key) #change the float key to a string
    #new_key = re.sub('\..*$','',new_key) #remove the decimal and everything following
    source_dict_timekeys[new_key] = source_dict_timekeys.pop(key) #new key on left of equal sign, old key on right

source_timekeys = ColumnDataSource(data=source_dict_timekeys)
plot2 = Figure(title='ion counts vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=950,plot_height=300)
plot2.vbar(x='x', bottom=0, width=0.5, top='y',color='firebrick',source=source_timekeys)

callback7 = CustomJS(args=dict(source=source_timekeys), code="""
     var data = source.data;

     //The object from the event ('doubletap') is stored as a cb_obj.
     //   The attribute x of cb_obj (the x coordinate double-clicked) is of interest
     var f = cb_obj.x

     //The retention times are the keys of the map object, data
     ris = Object.keys(data)

     //The keys are strings (as necessary) so they are converted to float numbers here for later mathamatical calculations
     var ris2 = []
     for (i = 0; i < ris.length; i++) {
         ris2[i] = parseFloat(ris[i])
     }

     //Some of the keys are not retention times and are variable names, they are NaN when converted to floats
     //   Those are removed here
     var ris2_real = []
     for (i = 0; i < ris2.length; i++) {
         if (ris2[i]==ris2[i]) {ris2_real.push(ris2[i])}
     }

     //An array with the absolute value of the difference between each retention index and the clicked x value is calculated
     var diff = ris2_real.map( function(value) {return value-f})
     var abs = diff.map( function(value) {return Math.abs(value)})

     //The minimum value of these differences is found
     //    The ... operator is needed to find the minimum of an array for some reason
     var min = Math.min(...abs)

     //The index of that minimum is found
     var ind = abs.indexOf(min)

     //The key of the retention time that should be plotted is found by finding the key closest to the double-click
     //    It must be a string
     new_key = ris2_real[ind]
     new_key = new_key.toString()

     //allows for debugging javascript
     //document.write(new_key)

     //set the x data
     x = data['x']
     y = data['y']

     //update to new y data
     to_change_to = data[new_key]
     for (i = 0; i < x.length; i++) {
         y[i] = to_change_to[i]
     }
     source.change.emit();
 """)

#indicate to change the second plot on a double-click on the first plot
plot.js_on_event('doubletap', callback7)

# Set up layouts and add to document
text_box0 = widgetbox(mz_text[0], width=175, height=20)
text_box1 = widgetbox(mz_text[1], width=175, height=20)
text_box2 = widgetbox(mz_text[2], width=175, height=20)
text_box3 = widgetbox(mz_text[3], width=175, height=20)
#layout = row(text_boxes, plot, plot2)

l = layout([
  [text_box0,text_box1,text_box2,text_box3],
  [plot],
  [plot2],
], sizing_mode='fixed')

#set the name shown in the tab to the sample_name
show(l)

#reset the bokeh plot ColumnDataSource object because it stores things dependent on the python instance
#    thus if you run this script multiple times in the same instance of python, the produced html file will grow linearly in size
#    every time you run it, it stores another copy of the data in the html file, unless the output is reset as below
reset_output()
