import numpy as np
import pickle #allows for saving of dictionary files
import copy

from bokeh.layouts import column, row, widgetbox
from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput
from bokeh.plotting import Figure, output_file, show

mz_plot = ['mz174','mz175','mz176','mz177','mz178']
mz_colors = ['red','blue','green','purple','orange']
file_directory = '/Users/nate/Desktop/netcdf_test/'
filename = 'tbdms01_t47d_wt.CDF'
output_data_file = file_directory + 'processed_data.p'
output_plot_file = file_directory + 'plot1.html'

output_file(output_plot_file)
plot=Figure(title='ion counts vs. time', x_axis_label='time (s)',y_axis_label='ion counts',plot_width=800,plot_height=400)

file_object = open(output_data_file,'rb')
file_data = pickle.load(file_object)
file_object.close()

source = {}
mz_text = {}
legend = {}
x_data = file_data[filename]['sats'] #this does not change with mz so it is set outside the loop

source_dict = file_data[filename]['ics_smooth_bc']
source_dict_keys = list(source_dict.keys())
for key in source_dict_keys:
    source_dict['mz'+str(int(key))] = source_dict.pop(key)

source_dict['x'] = x_data
source_dict['y0'] = source_dict[mz_plot[0]] #the ics for the current mz
source_dict['y1'] = source_dict[mz_plot[1]]
source_dict['y2'] = source_dict[mz_plot[2]]
source_dict['y3'] = source_dict[mz_plot[3]]
source_dict['y4'] = source_dict[mz_plot[4]]
source = ColumnDataSource(data=source_dict) #for bokeh widgets it is stored in a ColmnDataSource object

y = ['y0','y1','y2','y3','y4']

for j in [0,1,2,3,4]:
    plot.line('x',y[j],source=source,color=mz_colors[j]) #update the plot object for the current mz
    mz_text[j] = TextInput(title=mz_colors[j], value=str(mz_plot[j])) #the textbox widget, the value must be a string


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

callback4 = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    x = data['x']
    y4 = data['y4']
    to_change_to = data[f]
    for (i = 0; i < x.length; i++) {
        y4[i] = to_change_to[i]
    }
    source.change.emit();
""")

#inform bokeh to update the plot based on textbox values
mz_text[0].js_on_change('value', callback0)
mz_text[1].js_on_change('value', callback1)
mz_text[2].js_on_change('value', callback2)
mz_text[3].js_on_change('value', callback3)
mz_text[4].js_on_change('value', callback4)

# Set up layouts and add to document
text_boxes = widgetbox(mz_text[0],mz_text[1],mz_text[2],mz_text[3],mz_text[4],width=200,height=400)
layout = row(text_boxes, plot)

show(layout)
