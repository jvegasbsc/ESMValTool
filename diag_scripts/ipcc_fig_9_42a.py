"""
#############################################################################
ipcc_fig_9_42a.py
Author: Manuel Schlund (DLR, Germany)
ESMVal project
#############################################################################

Description
    Calculate and plot the equilibrium climate sensitivity (ECS) vs. the
    global mean surface temperature (GMSAT) for several CMIP5 models
    (see IPCC AR5 WG1 ch. 9, fig. 9.42a)

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
    20171107-A_schl_ma: written

#############################################################################
"""


# Basic python packages
import sys

# NetCDF4
import netCDF4 as nc

# ESMValTool python packages
from esmval_lib import ESMValProject


def main(project_info):
    """
    This is the main routine

    Parameters
    ----------
    project_info : dict
        Dictionary with project information
    """

    # Create instance of ESMValProject wrapper
    E = ESMValProject(project_info)

    config_file = E.get_configfile()    # Config files
    plot_dir = E.get_plot_dir()         # Plot directory
    datakeys = E.get_currVars()         # Current variables
    verbosity = E.get_verbosity()       # Verbosity state

    # Get model information for desired variable 'tas'
    tas_models = {}
    for datakey in datakeys:
        if (datakey == "tas"):
            tas_key = datakey
            tas_models = E.get_clim_model_filenames(tas_key)
    print(tas_models)                   # Print for dev. purposes


    # Iterate over all models
    for model in tas_models:
        pass

        # Initiallize model specific file
        tas_file = nc.Dataset(tas_models[model], 'r')

        # Get units
        tas_units = tas_file.variables[tas_key].units
        print(tas_units)
