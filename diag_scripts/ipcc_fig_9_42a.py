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
from auxiliary import error, info
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

    # Create instance of ESMValProject wrapper
    E = ESMValProject(project_info)

    # Get important information
    diag_name = E.get_diag_script_name()
    config_file = E.get_configfile()    # Config files
    plot_dir = E.get_plot_dir()         # Plot directory
    vars = E.get_currVars()             # Current variables
    verbosity = E.get_verbosity()       # Verbosity state

    # Check if all needed varibles are present
    if ("tas-degC" not in vars):
        error("no data for 'tas-degC' available, please check your namelist")
    if ("tas" not in vars):
        error("no data for 'tas' available, please check your namelist")
    if ("rtmt" not in vars):
        error("no data for 'rtmt' available, please check your namelist")

    # Write references:
    E.write_references(diag_name,       # diagnostic script name
                       ["A_schl_ma"],   # authors
                       [""],            # contributors
                       ["D_andrews12"], # diagnostic
                       [""],            # observations
                       ["P_esmval"],    # project
                       project_info,
                       verbosity,
                       False)

    # Print for dev. purposes
    print("config file: {0}".format(config_file))
    print("plot directory: {0}".format(plot_dir))
    print("variables: {0}".format(vars))
    print("verbosity: {0}".format(verbosity))
    print("")

    # Read configuration file
    modelconfig = ConfigParser.ConfigParser()
    modelconfig.read(config_file)
    area = modelconfig.get("test", "area")

    # Get all models
    models = E.get_all_clim_models()


    ###########################################################################
    # ECS calculation
    ###########################################################################

    # Dictionaries which will collect the data
    tas_piC = {}
    tas_4xCO2 = {}
    tas_delta = {}
    rtmt_piC = {}
    rtmt_4xCO2 = {}
    rtmt_delta = {}

    # Iterate over all models
    for model_path in models:
        model_info = models[model_path]

        # Try to get needed model information
        model_var = None
        model_name = None
        model_exp = None
        try:
            model_var = model_info["var"]
            model_name = model_info["name"]
            model_exp = model_info["experiment"]
        except KeyError:
            info("Could not retrieve all desired model information of " +
                 "model {0}".format(model_info), verbosity, 0)
            continue
        except:
            error("Unknown error while retrieving desired model information")

        # tas
        if (model_var == "tas"):
            GO = GridOperations(model_path, "tas")
            avg1 = GO._spatial_average(region=[[0,20],[70,90]])
            avg2 = GO._temporal_average(period="monthly")
            # print(avg1)
            # print(avg2)

        """
        print("---------MODEL: {0}-----------".format(model_name + "_" + \
                                                     model_exp + "_" + \
                                                     model_var))

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
        """

    print("\n****************************************************************")
