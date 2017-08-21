def ImportModules():

    import pdb #python debugger
    import os #allows to get working directory
    import re #allows regular expression use
    import plotly #for plotting, not used here?
    import importlib #allows fresh importing of modules
    import peakutils #allows for peak and baseline finding
    import bokeh.plotting as bkp #allows for making interactive plots with bokeh
    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas #a module which allows for making data frames
    import copy
    import scipy
    import savitzky_golay
    importlib.reload(savitzky_golay) #reload the custom function in case it was changed
    import FindBorders
    importlib.reload(FindBorders)
    import ExtendBounds
    importlib.reload(ExtendBounds)
    import fragment_library #a custom function
    importlib.reload(fragment_library) #reload the custom function in case it was changed
    import organize_ms_data
    importlib.reload(organize_ms_data)

    return(pdb,os,re,importlib,plotly,importlib,peakutils,bkp,np,pandas,copy,scipy,savitzky_golay,FindBorders,ExtendBounds,fragment_library,organize_ms_data)
