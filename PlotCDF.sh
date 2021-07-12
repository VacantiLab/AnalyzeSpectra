#!/bin/bash

source activate gcms10
export PYTHONPATH="${PYTHONPATH}:/Users/nate/git_hub_projects"

bokeh serve --show PlotCDF.py --port $1
