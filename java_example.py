def java_example():

    from bokeh.layouts import column
    from bokeh.models import CustomJS, ColumnDataSource, Slider, TextInput
    from bokeh.plotting import Figure, output_file, show
    import math

    output_file("java_example.html")

    x = [x*0.005 for x in range(0, 200)]
    to_plot = x
    z = [math.sin(15*item) for item in x]
    q = [math.sin(20*item) for item in x]

    source = ColumnDataSource(data=dict(x=x, to_plot=to_plot, z=z, q=q))

    plot = Figure(plot_width=400, plot_height=400)
    plot.line('x', 'to_plot', source=source, line_width=3, line_alpha=0.6)


    #somehow everytime this is run in the same python instance, the code below is added to the html file
    #it appears to be something CustomJS does to the python environment because it is stopped by deleting the html file
    #only starting a new instance of python stops the repitition of the below code
    callback = CustomJS(args=dict(source=source), code="""
        var data = source.data;
        var f = cb_obj.value
        x = data['x']
        to_plot = data['to_plot']
        to_change_to = data[f]
        for (i = 0; i < x.length; i++) {
            to_plot[i] = to_change_to[i]
        }
        source.change.emit();
    """)

    #The input to the text box is the name of the dictionary key containing the data you want to change to!
    text = TextInput(title='input', value='none')
    text.js_on_change('value', callback)

    layout = column(text, plot)

    show(layout)

    return()
