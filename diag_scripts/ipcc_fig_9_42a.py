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
        plot : Switch to plot the linear regression needed for the ECS
               calculation

Optional diag_script_info attributes (diagnostic specific)
    [main_plot]
        fontsize : Fonzsize used in the plot
    [ecs_plots]
        fontsize : Fontsize used in the plot
        xmin     : Left boundary of the plot
        xmax     : Right boundary of the plot
        ymin     : Lower boundary of the plot
        ymax     : Upper boundary of the plot

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
from auxiliary import info, warning, error
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
    global_conf = E.get_global_conf()
    diag_name = E.get_diag_script_name()
    vars = E.get_currVars()
    config_file = E.get_configfile()
    work_dir = E.get_work_dir()
    plot_dir = E.get_plot_dir()
    plot_file_type = E.get_graphic_format()
    write_plots = E.get_write_plots()
    verbosity = E.get_verbosity()
    exit_on_warning = E.get_exit_on_warning()
    if (plot_file_type not in plt.gcf().canvas.get_supported_filetypes()):
        info("Warning: selected file type for plots is not supported",
             verbosity, 0)
        plot_file_type = "ps"

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
    models = E.get_all_clim_models([var for var in VARIABLES])


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

        # Skip unnecesseary calculations
        if (model_var not in VARIABLES):
            continue
        if (model_exp not in VARIABLES[model_var]):
            continue
        grid_op = GridOperations(model_path, model_var, global_conf)
        info("Retrieving '{0}' from model '{1}' [{2}]".format(model_var,
                                                              model_name,
                                                              model_exp),
             verbosity, 1)
        units.update({model_var: grid_op.get_var_units()})

        # tas-degC
        if (model_var==TASDEGC):
            gmsat = grid_op.average(spatial_axis="all", period="total",
                                    spatial_weighting=True, region="global")[0]
            tasdegC_data[model_exp].update({model_name: gmsat})

        # tas
        if (model_var==TAS):
            tas = grid_op.average(spatial_axis="all", period="annual",
                                  spatial_weighting=True, region="global")
            tas_data[model_exp].update({model_name: tas})

        # rtmt
        if (model_var==RTMT):
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
    E.compare_models(tas_data[PICONTROL], tas_data[ABRUPT4XCO2], TAS+"/"+RTMT)
    E.compare_models(tasdegC_data[PICONTROL], tasdegC_data[HISTORICAL],
                     TASDEGC)

    # GMSAT
    gmsat_piC = tasdegC_data[PICONTROL]
    gmsat_hist = tasdegC_data[HISTORICAL]

    # ECS calculation (cf. Andrews et al. 2015)
    tas_piC = tas_data[PICONTROL]
    tas_4xCO2 = tas_data[ABRUPT4XCO2]
    rtmt_piC = rtmt_data[PICONTROL]
    rtmt_4xCO2 = rtmt_data[ABRUPT4XCO2]
    ecs_data = {}

    # Matplotlib instance
    fig, axes = plt.subplots()

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
            if (not write_plots):
                continue
            ecs_section = "ecs_plots"
            create_plot = E.get_config_option(modelconfig, ecs_section,
                                              "plot", False)
            if (create_plot is True):
                ecs_dir = plot_dir + "ecs_reg/"
                E.ensure_directory(ecs_dir)
                info("Create ECS regression plot for " + \
                     "model '{0}' in '{1}'".format(model, ecs_dir),
                     verbosity, 1)

                # Get values from configuration file
                cfg_options = {"fontsize": 16.0,
                               "xmin": 0.0, "xmax": 7.0,
                               "ymin": -2.0, "ymax": 10.0}
                cfg = E.get_config_options(modelconfig, ecs_section,
                                           cfg_options)

                # Get x and y of regression line
                x_reg = [cfg["xmin"]-1, cfg["xmax"]+1]
                y_reg = [reg_stats.slope*(cfg["xmin"]-1) + reg_stats.intercept,
                         reg_stats.slope*(cfg["xmax"]+1) + reg_stats.intercept]

                # Plot
                axes.plot(delta_tas, delta_rtmt, linestyle="none",
                          markeredgecolor="blue", markerfacecolor="none",
                          marker="s")
                axes.plot(x_reg, y_reg, color="black", linestyle="-")

                # Options
                axes.set_title(model, size=cfg["fontsize"]+2)
                axes.set_xlabel(TAS + " / " + units[TAS],
                                size=cfg["fontsize"])
                axes.set_ylabel(RTMT + " / " + units[RTMT],
                                size=cfg["fontsize"])
                axes.set_xlim(cfg["xmin"], cfg["xmax"])
                axes.set_ylim(cfg["ymin"], cfg["ymax"])
                axes.tick_params(labelsize=cfg["fontsize"]-2)
                axes.axhline(linestyle="dotted", c="black")
                axes.text(cfg["xmin"]+0.1, cfg["ymin"]+0.5,
                          "r = {:.2f},  ".format(reg_stats.rvalue) + \
                          r"$\alpha$ = {:.2f},  ".format(-reg_stats.slope) + \
                          "F = {:.2f},  ".format(reg_stats.intercept) + \
                          "ECS = {:.2f}".format(-reg_stats.intercept / \
                                                (2*reg_stats.slope)),
                          size=cfg["fontsize"])

                # Save plot
                filename = model + "." + plot_file_type
                fig.tight_layout()
                fig.savefig(ecs_dir + filename, orientation="landscape")
                axes.cla()

    # Empty line
    info("", verbosity, 1)


    ###########################################################################
    # Plot data
    ###########################################################################

    E.compare_models(ecs_data, gmsat_piC, "ECS/GMSAT[piControl]")
    E.compare_models(ecs_data, gmsat_hist, "ECS/GMSAT[historical]")
    piC_data = {}
    hist_data = {}

    # Collect data
    for model in ecs_data:
        if (model in gmsat_piC):
            piC_data.update({model: [ecs_data[model], gmsat_piC[model]]})
        if (model in gmsat_hist):
            hist_data.update({model: [ecs_data[model], gmsat_hist[model]]})

    # Plot data
    if (write_plots):
        info("Create IPCC AR5 WG1 fig. 9.42a plot in {0}".format(plot_dir),
             verbosity, 1)

        # Get values from configuration file
        main_plot_section = "main_plot"
        cfg_options = {"fontsize": 16.0}
        cfg = E.get_config_options(modelconfig, main_plot_section, cfg_options)

        # piControl
        for model in piC_data:
            style = E.get_model_style(model)
            axes.plot(piC_data[model][0], piC_data[model][1], linestyle="none",
                      markeredgecolor=style["color"],
                      markerfacecolor=style["facecolor"], marker=style["mark"],
                      markersize=cfg["fontsize"]-8, label="_"+model)

        # historical
        for model in hist_data:
            style = E.get_model_style(model)
            axes.plot(hist_data[model][0], hist_data[model][1],
                      linestyle="none", markeredgecolor=style["color"],
                      markerfacecolor=style["facecolor"], marker=style["mark"],
                      markersize=cfg["fontsize"]-4, label=model)

        # Options
        axes.set_title("IPCC AR5 WG1 - Fig. 9.42a", size=cfg["fontsize"]+2)
        axes.set_xlabel(r"ECS / $^\circ$C", size=cfg["fontsize"])
        axes.set_ylabel(r"GMSAT / $^\circ$C", size=cfg["fontsize"])
        axes.set_xlim(1.5, 5.0)
        axes.tick_params(labelsize=cfg["fontsize"]-2)
        legend = axes.legend(loc="upper left", fontsize=cfg["fontsize"],
                             bbox_to_anchor=(1.05, 1.0), borderaxespad=0.0)

        # Save plot
        filename = "IPCC-AR5-WG1_fig9-42a." + plot_file_type
        fig.savefig(plot_dir + filename, additional_artists=[legend],
                    bbox_inches="tight", orientation="landscape")
        plt.close("all")

    # Empty line
    info("", verbosity, 1)

    # Print data
    info("RESULTS:", verbosity, 0)
    for model in ecs_data:
        try:
            info("{0}:\tECS = {1}\tT[piC] = {2}\tT[hist] = {3}".format(
                model, ecs_data[model], gmsat_piC[model], gmsat_hist[model]),
                verbosity, 0)
        except KeyError:
            info("{0}:\tECS = {1}\tT[piC] = n.a.\tT[hist] = n.a.".format(
                model, ecs_data[model]), verbosity, 0)
