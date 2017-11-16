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
    for key in project_info:
        print("{0}: {1}\n".format(key, project_info[key]))

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

    # Get models
    models = E.get_cmip_clim_models()

    # Iterate over all models
    for model in models:
        model_name = models[model][0]["name"]
        model_var = models[model][0]["var"]
        model_path = models[model][1]

        print("--------------MODEL: {0}------------------".format(model_name))

        # Initiallize model specific file
        file = nc.Dataset(model_path, 'r')

        # Get data (tas)
        variables = file.variables
        data = variables[model_var][:]
        lat = variables["lat"][:]
        lon = variables["lon"][:]
        time = variables["time"][:]
        unit = file.variables[model_var].units

        # Get mean surface temp of every month
        GridOps = GridOperations(data, time, lat, lon)
        avg = GridOps.spatial_average()
        print("**** TOTAL AVERAGE: {0} = {1} {2} ****\n".format(model_var,
                                                              np.mean(avg),
                                                              unit))

        # Get plotting information
        color, dashes, width = E.get_model_plot_style(model_name)

    print("\n****************************************************************")
