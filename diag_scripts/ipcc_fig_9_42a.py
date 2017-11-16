# -*- coding: utf-8 -*-

"""
###############################################################################
ipcc_fig_9_42a.py
Author: Manuel Schlund (DLR, Germany)
ESMVal project
###############################################################################

Description
    Calculate and plot the equilibrium climate sensitivity (ECS) vs. the global
    mean surface temperature (GMSAT) for several CMIP5 models (see IPCC AR5 WG1
    ch. 9, fig. 9.42a).

Required diag_script_info attributes (diagnostics specific)
    none

Optional diag_script_info attributes (diagnostic specific)
    none

Required variable_info attributes (variable specific)
    none

Required variable attributes (defined in namelist)
    none

Caveats

Modification history
    20171109-A_schl_ma: written

###############################################################################
"""


# ESMValTool python packages
from esmval_lib import ESMValProject
from grid_operations import GridOperations

# NetCDF4
import netCDF4 as nc

# Basic python packages
import ConfigParser
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys


def main(project_info):
    """
    Arguments
        project_info : Dictionary containing project information

    Description
        This is the main routine of the diagnostic.
    """

    # Highlight user output
    print("\n****************************************************************")

    # Print project info
    print("PROJECT: {0}\n".format(project_info))

    # Create instance of ESMValProject wrapper
    E = ESMValProject(project_info)

    # Get important information
    config_file = E.get_configfile()    # Config files
    plot_dir    = E.get_plot_dir()      # Plot directory
    datakeys    = E.get_currVars()      # Current variables
    verbosity   = E.get_verbosity()     # Verbosity state
    experiment  = "equatorial"          # Current experiment

    # Print for dev. purposes
    print("config file: {0}".format(config_file))
    print("plot directory: {0}".format(plot_dir))
    print("datakeys: {0}".format(datakeys))
    print("verbosity: {0}".format(verbosity))

    # Read configuration file
    modelconfig = ConfigParser.ConfigParser()
    modelconfig.read(config_file)
    area = modelconfig.get(experiment, "area")

    # Get model information for desired variable 'tas'
    tas_models = {}
    for datakey in datakeys:
        if (datakey == "tas"):
            tas_key = datakey
            tas_models = E.get_clim_model_filenames(tas_key)
    print("models: {0}\n".format(tas_models)) # Print for dev. purposes

    # Iterate over all models
    for model in tas_models:

        # Initiallize model specific file
        tas_file = nc.Dataset(tas_models[model], 'r')

        # Get variables of file
        tas_variables = tas_file.variables
        print("VARIABLES OF FILE: {0}\n".format(tas_variables))
        print("DIMENSIONS OF FILE: {0}\n".format(tas_file.dimensions))

        data = tas_variables["tas"][:]
        print("DATA: {0}".format(data))
        print("SHAPE: {0}".format(data.shape))
        print("TYPE: {0}".format(type(data)))
        print("DIM: {0}\n".format(np.ndim(data)))

        lat = tas_variables["lat"][:]
        print("LAT: {0}\n".format(lat))
        lon = tas_variables["lon"][:]
        print("LON: {0}\n".format(lon))
        time = tas_variables["time"][:]
        print("TIME: {0}\n".format(time))

        # Get mean surface temp of every month
        GridOps = GridOperations(data, time, lat, lon)
        print("TOTAL MEAN = {0}".format(GridOps.spatial_average()))
        print("LAT MEAN = {0}".format(GridOps.spatial_average(axis="lat")))
        print("LON MEAN = {0}".format(GridOps.spatial_average(axis="lon")))

        # Get units
        tas_units = tas_file.variables[tas_key].units
        print("units: {0}".format(tas_units))

        # Get plotting information
        color, dashes, width = E.get_model_plot_style(model)
        print("color: {0}".format(color))
        print("dashes: {0}".format(dashes))
        print("width: {0}\n".format(width))

    print("\n****************************************************************")
