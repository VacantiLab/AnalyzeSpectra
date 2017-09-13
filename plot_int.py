import numpy as np
import pickle #allows for saving of dictionary files
import copy
import pdb

from bokeh.layouts import column, row, widgetbox, layout
from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput, Select
from bokeh.plotting import Figure, output_file, show, reset_output
from bokeh.io import curdoc

mz_plot = ['174.0','175.0','176.0','177.0']
time_key_plot = 't738.54700000000003'
mz_colors = ['red','blue','green','purple']
file_directory = '/Users/nate/Desktop/netcdf_test/'
filename = 'tbdms01_t47d_wt.CDF'
output_data_file = file_directory + 'processed_data.p'
output_plot_file = file_directory + 'plot1.html'

output_file(output_plot_file, title=filename) #title sets the title shown on the browser tab
plot=Figure(title='ion counts vs. time', x_axis_label='time (s)',y_axis_label='ion counts',plot_width=950,plot_height=350)

file_object = open(output_data_file,'rb')
file_data = pickle.load(file_object)
file_object.close()

source = {}
mz_text = {}
legend = {}

x_data = file_data[filename]['sats'] #this does not change with mz so it is set outside the loop
blank_data = np.zeros(len(x_data))
for i in range(0,len(blank_data)):
    blank_data[i] = np.nan


#change the key names to strings
source_dict = file_data[filename]['ics_smooth_bc']
source_dict_keys = list(source_dict.keys())
for key in source_dict_keys:
    source_dict[str(key)] = source_dict.pop(key) #new key on left of equal sign, old key on right

#make smaller for testing
#new_source_dict = {}
#new_source_dict['174.05'] = source_dict['174.05']
#new_source_dict['175.05'] = source_dict['175.05']
#new_source_dict['176.05'] = source_dict['176.05']
#new_source_dict['177.05'] = source_dict['177.05']
#new_source_dict['178.05'] = source_dict['178.05']
#source_dict = copy.copy(new_source_dict)

source_dict['x'] = x_data
source_dict['blank'] = blank_data
source_dict['y0'] = source_dict[mz_plot[0]] #the ics for the current mz
source_dict['y1'] = source_dict[mz_plot[1]]
source_dict['y2'] = source_dict[mz_plot[2]]
source_dict['y3'] = source_dict[mz_plot[3]]
source = ColumnDataSource(data=source_dict) #for bokeh widgets it is stored in a ColmnDataSource object

y = ['y0','y1','y2','y3']

for j in [0,1,2,3]:
    plot.line('x',y[j],source=source,color=mz_colors[j]) #update the plot object for the current mz
    #mz_text[j] = TextInput(title=mz_colors[j], value=str(mz_plot[j])) #the textbox widget, the value must be a string
    mz_text[j] = Select(title=mz_colors[j], value=str(mz_plot[j]), options=list(source_dict.keys()))

callback0 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y0 = data['y0']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y0[i] = to_change_to[i]
    }
    source.change.emit();
""")

callback1 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y1 = data['y1']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y1[i] = to_change_to[i]
    }
    source.change.emit();
""")

callback2 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y2 = data['y2']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y2[i] = to_change_to[i]
    }
    source.change.emit();
""")

callback3 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y3 = data['y3']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y3[i] = to_change_to[i]
    }
    source.change.emit();
""")

#inform bokeh to update the plot based on textbox values
mz_text[0].js_on_change('value', callback0)
mz_text[1].js_on_change('value', callback1)
mz_text[2].js_on_change('value', callback2)
mz_text[3].js_on_change('value', callback3)

#Make the ion-count vs. mz plot for each scan acquisition time
source_dict_timekeys = file_data[filename]['ics_smooth_timekeys']
source_dict_timekeys_keys = list(source_dict_timekeys.keys())
for key in source_dict_timekeys_keys:
    source_dict_timekeys[str(key)] = source_dict_timekeys.pop(key)

#make smaller for testing
#source_dict_timekeys_new = {}
#source_dict_timekeys_new[list(source_dict_timekeys.keys())[0]] = source_dict_timekeys[list(source_dict_timekeys.keys())[0]]
#source_dict_timekeys_new[list(source_dict_timekeys.keys())[1]] = source_dict_timekeys[list(source_dict_timekeys.keys())[1]]
#source_dict_timekeys_new[list(source_dict_timekeys.keys())[2]] = source_dict_timekeys[list(source_dict_timekeys.keys())[2]]
#source_dict_timekeys = copy.copy(source_dict_timekeys_new)

source_dict_timekeys['x'] = file_data[filename]['mz_vals']
test_time_value = list(source_dict_timekeys.keys())[0]
source_dict_timekeys['y'] = source_dict_timekeys[test_time_value]

source_timekeys = ColumnDataSource(data=source_dict_timekeys)
plot2 = Figure(title='ion counts vs. mz', x_axis_label='m/z',y_axis_label='ion counts',plot_width=950,plot_height=350)
plot2.vbar(x='x', bottom=0, width=0.5, top='y',color='firebrick',source=source_timekeys)

callback6 = CustomJS(args=dict(source=source_timekeys), code="""
     var data = source.data;
     var f = cb_obj.value
     f_s = f.toString()
     x = data['x']
     y = data['y']
     to_change_to = data[f_s]
     for (i = 0; i < x.length; i++) {
         y[i] = to_change_to[i]
     }
     source.change.emit();
 """)

time_select = Select(title="Scan Acquisition Time:", value=test_time_value, options=list(source_dict_timekeys.keys()))
time_slider = Slider(start=395.572, end=2716.746, value=395.572, step=0.343, title='Scan Acquisition Time', callback=callback6, callback_policy='mouseup')

callback5 = CustomJS(args=dict(source=source_timekeys), code="""
     var data = source.data;
     var f = cb_obj.value
     x = data['x']
     y = data['y']
     to_change_to = data[f]
     for (i = 0; i < x.length; i++) {
         y[i] = to_change_to[i]
     }
     source.change.emit();
 """)

time_select.js_on_change('value', callback5)
#time_slider.js_on_change('value', callback6)

# Set up layouts and add to document
text_box0 = widgetbox(mz_text[0], width=100, height=20)
text_box1 = widgetbox(mz_text[1], width=100, height=20)
text_box2 = widgetbox(mz_text[2], width=100, height=20)
text_box3 = widgetbox(mz_text[3], width=100, height=20)
time_box = widgetbox(time_select,width=200,height=350)
time_slider_box = widgetbox(time_slider,width=1000,height=30)
#layout = row(text_boxes, plot, plot2)

l = layout([
  [text_box0,text_box1,text_box2,text_box3],
  [plot],
  [time_slider_box],
  [plot2],
], sizing_mode='fixed')

#set the name shown in the tab to the filename
show(l)

#reset the bokeh plot ColumnDataSource object because it stores things dependent on the python instance
#    thus if you run this script multiple times in the same instance of python, the produced html file will grow linearly in size
#    every time you run it, it stores another copy of the data in the html file, unless the output is reset as below
reset_output()
