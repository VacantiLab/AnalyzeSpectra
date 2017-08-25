import numpy as np
import pickle #allows for saving of dictionary files
import copy

from bokeh.layouts import column, row, widgetbox, layout
from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput, Select
from bokeh.plotting import Figure, output_file, show, reset_output

mz_plot = ['174','175','176','177','178']
time_key_plot = 't738.54700000000003'
mz_colors = ['red','blue','green','purple','orange']
file_directory = '/Users/nate/Desktop/netcdf_test/'
filename = 'tbdms01_t47d_wt.CDF'
output_data_file = file_directory + 'processed_data.p'
output_plot_file = file_directory + 'plot1.html'

output_file(output_plot_file)
plot=Figure(title='ion counts vs. time', x_axis_label='time (s)',y_axis_label='ion counts',plot_width=800,plot_height=350)

file_object = open(output_data_file,'rb')
file_data = pickle.load(file_object)
file_object.close()

source = {}
mz_text = {}
legend = {}
x_data = file_data[filename]['sats'] #this does not change with mz so it is set outside the loop

#change the key names to strings
source_dict = file_data[filename]['ics_smooth_bc']
source_dict_keys = list(source_dict.keys())
for key in source_dict_keys:
    source_dict[str(int(key))] = source_dict.pop(key) #new key on left of equal sign, old key on right

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

#Make the ion-count vs. mz plot for each scan acquisition time
source_dict_timekeys = file_data[filename]['ics_smooth_timekeys']
source_dict_timekeys_keys = list(source_dict_timekeys.keys())
for key in source_dict_timekeys_keys:
    source_dict_timekeys[str(key)] = source_dict_timekeys.pop(key)

source_dict_timekeys['x'] = file_data[filename]['mz_vals']
test_time_value = list(source_dict_timekeys.keys())[0]
source_dict_timekeys['y'] = source_dict_timekeys[test_time_value]

source_timekeys = ColumnDataSource(data=source_dict_timekeys)
plot2 = Figure(title='ion counts vs. mz', x_axis_label='mz',y_axis_label='ion counts',plot_width=800,plot_height=350)
plot2.vbar(x='x', bottom=0, width=0.5, top='y',color='firebrick',source=source_timekeys)

time_select = Select(title="Scan Acquisition Time:", value=test_time_value, options=list(source_dict_timekeys.keys()))

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

# Set up layouts and add to document
text_boxes = widgetbox(mz_text[0],mz_text[1],mz_text[2],mz_text[3],mz_text[4],width=200,height=350)
time_box = widgetbox(time_select,width=200,height=350)
#layout = row(text_boxes, plot, plot2)

l = layout([
  [text_boxes, plot],
  [time_box,plot2],
], sizing_mode='fixed')

show(l)

#reset the bokeh plot ColumnDataSource object because it stores things dependent on the python instance
#    thus if you run this script multiple times in the same instance of python, the produced html file will grow linearly in size
#    every time you run it, it stores another copy of the data in the html file, unless the output is reset as below
reset_output()
