#/usr/bin/env python2
# -*- coding: utf-8 -*-


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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


class Plot2D(object):
    """
    Description
        Basic class for 2-dimensional plotting

    Contents
        method
    """

    LATS = ['latitude']
    LONS = ['longitude']
    TIME = ['time']
    SIDE = 5

    def __init__(self, cube, **kwargs):
        """
        Arguments
            cube : iris cube

            kwargs:
                x_line_plot : Add line plot below averaging over x-axis
                              ('top', 'bottom', None)
                y_line_plot : Add line plot on right averaging over y-axis
                              ('left', 'right', None)

                TODO:
                    swap x/y axis

        Description
            Initializes the class and performs the actual plotting.

        Modification history
            20180207-A_schl_ma: written
        """

        # Check arguments
        if not isinstance(cube, iris.cube.Cube):
            raise TypeError("Invalid input: expected iris cube")
        self.cube = iris.util.squeeze(cube)
        if (self.cube.ndim != 2):
            raise TypeError("Invalid input: expected 2-dimensional iris cube")

        # Process kwargs
        self.xline = None
        try:
            if (kwargs['x_line_plot'] is not None):
                self.xline = kwargs['x_line_plot']
        except KeyError:
            pass
        self.yline = None
        try:
            if (kwargs['y_line_plot'] is not None):
                self.yline = kwargs['y_line_plot']
        except KeyError:
            pass

        # Get dimension names
        dim_names = [dim.standard_name for dim in self.cube.dim_coords]
        for dim in dim_names:
            if (dim in self.__class__.LATS):
                self.lat_var = dim
                break
            else:
                self.lat_var = None
        for dim in dim_names:
            if (dim in self.__class__.LONS):
                self.lon_var = dim
                break
            else:
                self.lon_var = None
        for dim in dim_names:
            if (dim in self.__class__.TIME):
                self.time_var = dim
                break
            else:
                self.time_var = None

        # Lat/lon plot
        if (self.lat_var is not None and self.lon_var is not None):
            self.x = self.lon_var
            self.y = self.lat_var

        # Lat/time plot
        elif (self.lat_var is not None and self.time_var is not None):
            self.x = self.time_var
            self.y = self.lat_var

        # Lon/time plot
        elif (self.lat_var is not None and self.lon_var is not None):
            self.x = self.lon_var
            self.y = self.time_var

        # Default case
        else:
            raise TypeError("Invalid input: cube does not contain supported " +
                            "dimensions")

        self.__plot()

###############################################################################

    def __plot(self):
        """
        Arguments
            None

        Description
            Actual plotting routine

        Modification history
            20180207-A_schl_ma: written
        """

        # Contour plot
        G = gridspec.GridSpec(self.__class__.SIDE, self.__class__.SIDE)
        ax_main = plt.subplot(G[1:-1, 1:-1])
        qplt.contourf(self.cube)

        # If lat and lon, plot map
        if (self.x == self.lon_var and self.y == self.lat_var):
            plt.gca().coastlines()

        # Line plots
        if (self.xline == 'top'):
            ax_xline = plt.subplot(G[0, 1:-1])
        elif (self.xline == 'bottom'):
            ax_xline = plt.subplo(G[-1, 1:-1])
        elif (self.xline is None):
            ax_xline = None
        else:
            raise ValueError("Invalid input: x_line_plot needs to be " + \
                             "'top', 'bottom' or None")

        # Kwargs
        if False:
            if (kwargs['x_line_plot'] is not None):
                cube_line = cube.collapsed(lon, iris.analysis.MEAN)
                if (kwargs['x_line_plot'] == 'top'):
                    plt.subplot(332)

        plt.savefig(os.path.join(PLOT_DIR, PLOT_NAME))
        plt.clf()

Plot2D(latlon_cube, x_line_plot='top')
