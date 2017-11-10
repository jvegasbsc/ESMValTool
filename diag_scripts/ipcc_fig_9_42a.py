"""
###############################################################################
ipcc_fig_9_42a.py
Author: Manuel Schlund (DLR, Germany)
ESMVal project
###############################################################################

Description
    Calculate and plot the equilibrium climate sensitivity (ECS) vs. the global
    mean surface temperature (GMSAT) for several CMIP5 models (see IPCC AR5 WG1
    ch. 9, fig. 9.42a)

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
from latlon import gridcell_area

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
    This is the main routine of the diagnostic.

    Arguments
        project_info : Dictionary containing project information
    """

    # Highlight user output:
    print("\n****************************************************************")

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
        pass

        # Initiallize model specific file
        tas_file = nc.Dataset(tas_models[model], 'r')

        # Get variables of file
        tas_variables = tas_file.variables
        print("CONTENT OF FILE: {0}\n".format(tas_variables))
        # print("CONTENT OF TAS: {0}\n".format(tas_variables["tas"]))
        # print("TAS: {0}\n".format(tas_variables["tas"][:]))
        # print("TIME: {0}\n".format(tas_variables["time"][:]))

        # Get mean surface temp of every day


        # Get units
        tas_units = tas_file.variables[tas_key].units
        print("units: {0}".format(tas_units))

        # Get plotting information
        color, dashes, width = E.get_model_plot_style(model)
        print("color: {0}".format(color))
        print("dashes: {0}".format(dashes))
        print("width: {0}\n".format(width))

        # Get model data (equatorial experiment, Atlantic Ocean)
        tas_data = E.get_model_data(modelconfig,
                                    experiment,
                                    area,
                                    tas_key,
                                    tas_file)

        # Calculate annual means
        raw_data = tas_data[2]
        print("PURE DATA: {0}\n".format(raw_data))

        # Get average for every year
        annual_averages = E.average_data(raw_data, 0)
        print("ANNUAL AVERAGES: {0}\n".format(annual_averages))

        # Print model data
        print("SHAPE OF EXPERIMENT DATA: {0}".format(np.shape(tas_data)))
        print("LENGTH OF EXPERIMENT DATA: {0}\n".format(len(tas_data)))

        # Print example areas of grid:
        print("GRID(1,2,3,4) = {0}".format(gridcell_area(1,2,3,4)))
        print("GRID(90,92,45,42) = {0}".format(gridcell_area(90,92,45,42)))

    print("\n****************************************************************")
