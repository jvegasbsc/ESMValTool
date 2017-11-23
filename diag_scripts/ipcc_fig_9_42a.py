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

    ###########################################################################
    # Variables and experiments needed for this diagnostic
    ###########################################################################

    TASDEGC = "tas-degC"
    TAS = "tas"
    RTMT = "rtmt"
    HISTORICAL = "historical"
    PICONTROL = "piControl"
    ABRUPT4XCO2 = "abrupt4xCO2"

    VARIABLES = {TASDEGC: [HISTORICAL],
                 TAS: [PICONTROL, ABRUPT4XCO2],
                 RTMT: [PICONTROL, ABRUPT4XCO2]}


    ###########################################################################
    # Get namelist information
    ###########################################################################

    # Create instance of ESMValProject wrapper
    E = ESMValProject(project_info)

    # Get information
    diag_name = E.get_diag_script_name()
    config_file = E.get_configfile()
    plot_dir = E.get_plot_dir()
    vars = E.get_currVars()
    verbosity = E.get_verbosity()

    # Check if all needed varibles are present
    for var in VARIABLES:
        if (var not in vars):
            error("no data for variable '{0}' available, ".format(var) + \
                  "please check your namelist")

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

    # Read configuration file
    modelconfig = ConfigParser.ConfigParser()
    modelconfig.read(config_file)
    area = modelconfig.get("test", "area")

    # Get all models
    models = E.get_all_clim_models()


    ###########################################################################
    # Collect data of the models
    ###########################################################################

    # Dictionaries which will collect the data for all models
    tasdegC_data = {}
    tas_data = {}
    rtmt_data = {}

    # Iterate over all models
    for model_path in models:
        model_info = models[model_path]

        # Try to get necessary model information
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

        # Setup GridOperations member
        grid_op = GridOperations(model_path, model_var)

        # tas-degC:
        if (model_var==TASDEGC and model_exp in VARIABLES[TASDEGC]):
            print(model_name)
            print(model_var)
            print(model_exp)

            # Get GMSAT
            gmsat = grid_op.average(spatial_axis="all", period="total",
                                    spatial_weighting=True, region="global")
            print(gmsat)
            print("")

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
        avg = GridOps._spatial_average()
        print("**** TOTAL AVERAGE: {0} = {1} {2} ****\n".format(model_var,
                                                                np.mean(avg),
                                                                unit))

        # Get plotting information
        color, dashes, width = E.get_model_plot_style(model_name)
        """

    print("\n****************************************************************")
