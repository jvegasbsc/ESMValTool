#/usr/bin/env python2
# -*- coding: utf-8 -*-


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import iris
import iris.plot as iplt
import iris.quickplot as qplt


# Example dataset (simple 2D tas)
LATLON = '/home/schl_m8/ESMValTool/diag_scripts/aux/C3S_511/test_latlon.nc'
LATTIME = '/home/schl_m8/ESMValTool/diag_scripts/aux/C3S_511/test_lattime.nc'
LONTIME = '/home/schl_m8/ESMValTool/diag_scripts/aux/C3S_511/test_lontime.nc'
latlon_cube = iris.load_cube(LATLON)
lattime_cube = iris.load_cube(LATTIME)
lontime_cube = iris.load_cube(LONTIME)

# User arguments (TODO: substitute with class attributes):
PLOT_DIR = '/home/schl_m8/ESMValTool/diag_scripts/aux/C3S_511'
PLOT_NAME = 'latlon.pdf'


def plot_2D(cube, **kwargs):
    """
    2 dimensional plotting routine.
    """

    # User arguments (TODO: substitue with class attributes)
    plot_dir = '/home/schl_m8/ESMValTool/diag_scripts/aux/C3S_511'
    plot_name = 'latlon.pdf'

    # Function variables
    lat = 'lat'
    lon = 'lon'
    time = 'time'

    # Check arguments
    if not isinstance(cube, iris.cube.Cube):
        raise TypeError("Invalid input: expected iris cube")

    # Squeeze cube to remove dimensions of length 1
    cube = iris.util.squeeze(cube)

    # Get dimension names
    dim_names = [dim.var_name for dim in cube.dim_coords]

    # Setup matplotlib

    # Lat/lon plot
    if (lat in dim_names and lon in dim_names):
        __latlon_plot(cube)

    # Lat/time plot
    elif (lat in dim_names and time in dim_names):
        __lattime_plot(cube)

    # Lon/time plot
    elif (lon in dim_names and time in dim_names):
        __lontime_plot(cube)





    print(cube.dim_coords)


def __latlon_plot(cube):
    """
    Lat/lon plot.
    """
    qplt.contourf(cube)
    plt.gca().coastlines()
    plt.savefig(os.path.join(PLOT_DIR, PLOT_NAME))


def __lattime_plot(cube):
    """
    Lat/time plot.
    """
    pass


def __lontime_plot(cube):
    """
    Lon/time plot.
    """
    pass


plot_2D(latlon_cube)
