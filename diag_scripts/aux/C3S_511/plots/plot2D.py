#/usr/bin/env python2
# -*- coding: utf-8 -*-


import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cf_units as units


class Plot2D(object):
    """
    Description
        Basic class for 2-dimensional plotting

    Contents
        method
    """

    # Class attributes
    LATS = ['latitude']     # accepted lat names
    LONS = ['longitude']    # accepted lon names
    TIME = ['time']         # accepted time names
    GRIDSPEC_LENGTH = 5     # size of the matplotlib gridspec
    MPLSTYLE = 'diag_scripts/aux/C3S_511/plots/default.mplstyle'

    def __init__(self, cube, **kwargs):
        """
        Arguments
            cube : iris cube

            kwargs:
                summary_plot : Add summary line plot (True, False)

                TODO:
                    optional summary placement
                    check kwargs
                    swap x/y axis

        Description
            Initializes the class.

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
        if ('summary_plot' in kwargs):
            if (not isinstance(kwargs['summary_plot'], bool)):
                raise TypeError("Invalid input: summary_plot should be " +
                                "True or False")
            self.summary_plot = kwargs['summary_plot']
        else:
            self.summary_plot = False

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
            self.plot_type = 'latlon'
            if (self.summary_plot):
                self.summary_location = 'right'
                self.colorbar = 'horizontal'
            else:
                self.colorbar = 'vertical'

        # Lat/time plot
        elif (self.lat_var is not None and self.time_var is not None):
            self.plot_type = 'lattime'
            if (self.summary_plot):
                self.summary_location = 'right'
                self.colorbar = 'horizontal'
            else:
                self.colorbar = 'vertical'

        # Lon/time plot
        elif (self.lon_var is not None and self.time_var is not None):
            self.plot_type = 'lontime'
            if (self.summary_plot):
                self.summary_location = 'bottom'
            self.colorbar = 'vertical'

        # Default case
        else:
            raise TypeError("Invalid input: cube does not contain supported " +
                            "dimensions")

        # Setup matplotlib
        plt.style.use(self.__class__.MPLSTYLE)

###############################################################################

    def plot(self):
        """
        Arguments
            None

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180207-A_schl_ma: written
        """

        # Summary plot yes/no
        if (not self.summary_plot):
            G = gridspec.GridSpec(1, 1)
        elif (self.summary_location == 'right'):
            G = gridspec.GridSpec(2, 2, width_ratios=[4, 1],
                                  height_ratios=[5, 1])
        elif (self.summary_location == 'bottom'):
            G = gridspec.GridSpec(2, 2, width_ratios=[5, 1],
                                  height_ratios=[4, 1])

        # Main plot
        plt.subplot(G[0])
        if (self.plot_type == 'latlon'):
            x = self.cube.coord(self.lon_var).points
            y = self.cube.coord(self.lat_var).points
            z = self.cube.data
            collapse = self.lon_var
            line_plot_axes = 1
            cb_axes = 2
            m = bm.Basemap(projection='cyl', llcrnrlat=-90.0, urcrnrlat=90.0,
                           llcrnrlon=0.0, urcrnrlon=360.0, lat_ts=20.0,
                           resolution='c')
            m.drawcoastlines()
            x, y = m(*np.meshgrid(x, y))
            ax_main = m.contourf(x, y, z)
        else:
            if (self.plot_type == 'lattime'):
                x = self.cube.coord(self.time_var).points
                y = self.cube.coord(self.lat_var).points
                z = self.cube.data.T
                collapse = self.time_var
                line_plot_axes = 1
                cb_axes = 2
            elif (self.plot_type == 'lontime'):
                x = self.cube.coord(self.lon_var).points
                y = self.cube.coord(self.time_var).points
                z = self.cube.data
                collapse = self.time_var
                line_plot_axes = 2
                cb_axes = 1
            ax_main = plt.contourf(x, y, z)

        # Line plot
        if (self.summary_plot):
            ax_cb = plt.subplot(G[cb_axes])
            cb = plt.colorbar(ax_main, cax=ax_cb, orientation=self.colorbar)
            plt.subplot(G[line_plot_axes])
            if (collapse == self.lat_var):
                grid_areas = iris.analysis.cartography.area_weights(self.cube)
            else:
                grid_areas = None
            cube_line = self.cube.collapsed(collapse, iris.analysis.MEAN,
                                            weights=grid_areas)

            if (self.plot_type == 'latlon' or self.plot_type == 'lattime'):
                x = cube_line
                y = cube_line.coord(self.lat_var)
            else:
                x = cube_line.coord(self.lon_var)
                y = cube_line
            ax_line = iplt.plot(x, y)
        else:
            cb = plt.colorbar(ax_main, orientation=self.colorbar)

        plt.tight_layout()
        return plt.gcf()
