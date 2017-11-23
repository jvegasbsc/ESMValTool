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
    [ecs_plots]
        plot: Switch to plot the linear regression needed for the ECS calculation

Optional diag_script_info attributes (diagnostic specific)
    [ecs_plots]
        file_type: File type of the plots (default: 'eps')

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
from scipy import stats
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

    VARIABLES = {TASDEGC: [PICONTROL, HISTORICAL],
                 TAS: [PICONTROL, ABRUPT4XCO2],
                 RTMT: [PICONTROL, ABRUPT4XCO2]}


    ###########################################################################
    # Get namelist information
    ###########################################################################

    # Create instance of ESMValProject wrapper
    E = ESMValProject(project_info)

    # Get information
    diag_name = E.get_diag_script_name()
    vars = E.get_currVars()
    config_file = E.get_configfile()
    plot_dir = E.get_plot_dir()
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

    # Get all models
    models = E.get_all_clim_models()


    ###########################################################################
    # Collect data of the models
    ###########################################################################

    info("", verbosity, 1)
    info("Starting calculation", verbosity, 1)
    info("", verbosity, 1)

    # Dictionaries which will collect the data for all models
    tasdegC_data = {exp: {} for exp in VARIABLES[TASDEGC]}
    tas_data = {exp: {} for exp in VARIABLES[TAS]}
    rtmt_data = {exp: {} for exp in VARIABLES[RTMT]}
    units = {}

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

        # Setup GridOperations member
        grid_op = GridOperations(model_path, model_var)
        info("Retrieving '{0}' from model '{1}' [{2}]".format(model_var,
                                                              model_name,
                                                              model_exp),
             verbosity, 1)
        units.update({model_var: grid_op.get_var_units()})

        # tas-degC
        if (model_var==TASDEGC and model_exp in VARIABLES[TASDEGC]):
            gmsat = grid_op.average(spatial_axis="all", period="total",
                                    spatial_weighting=True, region="global")[0]
            tasdegC_data[model_exp].update({model_name: gmsat})

        # tas
        if (model_var==TAS and model_exp in VARIABLES[TAS]):
            tas = grid_op.average(spatial_axis="all", period="annual",
                                  spatial_weighting=True, region="global")
            tas_data[model_exp].update({model_name: tas})

        # rtmt
        if (model_var==RTMT and model_exp in VARIABLES[RTMT]):
            rtmt = grid_op.average(spatial_axis="all", period="annual",
                                   spatial_weighting=True, region="global")
            # Add to dictionary
            rtmt_data[model_exp].update({model_name: rtmt})

    # Empty line
    info("", verbosity, 1)


    ###########################################################################
    # Process data
    ###########################################################################

    # Check if data is missing
    missing_models = [model for model in tas_data[PICONTROL] if \
                      model not in tas_data[ABRUPT4XCO2]]
    for model in missing_models:
        info("Warning: model '{0}' does not contain data ".format(model) + \
             "of all needed experiments, check your namelist", verbosity, 0)

    # ECS calculation (cf. Andrews et al. 2015)
    tas_piC = tas_data[PICONTROL]
    tas_4xCO2 = tas_data[ABRUPT4XCO2]
    rtmt_piC = rtmt_data[PICONTROL]
    rtmt_4xCO2 = rtmt_data[ABRUPT4XCO2]
    ecs_data = {}

    # Iterate over all models (same for all variables)
    for model in tas_piC:
        if (model in tas_4xCO2):
            delta_tas = tas_4xCO2[model]-tas_piC[model]
            delta_rtmt = rtmt_4xCO2[model]-rtmt_piC[model]

            # Perform linear regression
            reg_stats = stats.linregress(delta_tas, delta_rtmt)
            ecs = -reg_stats.intercept / (2*reg_stats.slope)
            ecs_data.update({model: ecs})

            # Create plot for this linear regression if desired
            try:
                create_plot = modelconfig.getboolean("ecs_plots", "plot")
            except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                info("Warning: config file {0} does ".format(config_file) + \
                     "not contain option 'plot' in section 'ecs_plot'",
                     verbosity, 0)
                create_plot = False
            if (create_plot is True):
                ecs_dir = plot_dir + "ecs_reg/"
                E.ensure_directory(ecs_dir)
                info("Create ECS regression plot for " + \
                     "model '{0}'".format(model), verbosity, 1)

                # Plot
                fig = plt.figure()
                axes = fig.add_subplot(111)
                axes.scatter(delta_tas, delta_rtmt)
                axes.set_title(model)
                axes.set_xlabel(TAS + " / " + units[TAS])
                axes.set_ylabel(RTMT + " / " + units[RTMT])

                # Get file name
                default_filetype = "eps"
                try:
                    filetype = modelconfig.get("ecs_plots", "filetype")
                    if (filetype not in fig.canvas.get_supported_filetypes()):
                        filetype = default_filetype
                except (ConfigParser.NoSectionError,
                        ConfigParser.NoOptionError):
                    filetype = default_filetype
                filename = model + "." + filetype

                # Save plot
                fig.savefig(ecs_dir + filename)

    print(ecs_data)

    """
        # Get plotting information
        color, dashes, width = E.get_model_plot_style(model_name)
    """
