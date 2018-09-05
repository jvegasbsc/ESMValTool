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
import cf_units as unit
import matplotlib.cm as mpl_cm
import cartopy.crs as ccrs
#import random
import sys
from matplotlib.ticker import FuncFormatter
from string import ascii_lowercase

def label_in_perc(x, pos=0):
    return '%1.1f%%' % (x*100)

MPLSTYLE = os.path.dirname(os.path.realpath(__file__)) + os.sep + 'default.mplstyle'


class PlotHist(object):
    """
    Description
        Basic class for plotting histograms

    Contents
        method plot
    """

###############################################################################

    def __init__(self, data):
        """
        Arguments
            data : input data (iris cube or array)

        Description
            Initializes the class.

        Modification history
            20180209-A_schl_ma: written
        """

        # Check arguments
        if isinstance(data, iris.cube.Cube):
            try:
                self.name = data.long_name
            except:
                pass
            self.units = ' [' + str(data.units) + ']'
            self.data = np.ravel(data.data)
        elif isinstance(data, np.ndarray):
            self.name = 'data'
            self.units = ''
            self.data = data
        else:
            raise TypeError("Invalid input: expected iris cube or numpy array")

        # Matplotlib
        plt.style.use(MPLSTYLE)
        self.fig, self.ax = plt.subplots()

###############################################################################

    def plot(self, bins=20, x_label=None, y_label=None, title=None,
             color='green', alpha=0.75):
        """
        Arguments
            bins    : number of bins
            x_label : label of x-axis
            y_label : label of y-axis
            title   : title of the plot
            color   : color of the histogramm
            alpha   : transparency of the historgramm

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180209-A_schl_ma: written
        """

        # Parse arguments
        if (x_label is None):
            x_label = self.name + self.units
        if (y_label is None):
            y_label = 'frequency'

        # Create historgramm
#        self.ax.hist(self.data, bins, density=False, facecolor=color, alpha=alpha)
#        print self.ax.__dict__.keys()
        try:
            hist, bins = np.histogram(self.data.data[np.logical_not(self.data.mask)], bins=bins)
        except:
            hist, bins = np.histogram(self.data.data, bins=bins)
        hist = hist.astype(float)/hist.sum()
        binWidth = bins[1] - bins[0]
        self.ax.bar(bins[:-1]+0.5*binWidth, hist, binWidth, facecolor=color, alpha=alpha)
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        self.ax.yaxis.set_major_formatter(FuncFormatter(label_in_perc))
        if (title is not None):
            self.ax.set_title(title)

        self.fig.tight_layout()
        return self.fig


class PlotScatter(object):
    """
    Description
        Basic class for creating scatter plots

    Contents
        method plot
    """

###############################################################################

    def __init__(self, data1, data2):
        """
        Arguments
            data1 : single cube (y-axis)
            data2 : list of cubes which should be compared (x-axis)

        Description
            Initializes the class.

        Modification history
            20180209-A_schl_ma: written
        """

        # Check arguments
        if (not isinstance(data1, iris.cube.Cube)):
            raise TypeError("Invalid input: expected iris cube")
        if isinstance(data2, iris.cube.Cube):
            self.data2 = [data2]
        else:
            self.data2 = data2
        for d in self.data2:
            if (not isinstance(d, iris.cube.Cube)):
                raise TypeError("Invalid input: expected list of iris cubes " +
                                "or a single cube")
            if (data1.coords != d.coords):
                raise ValueError("Invalid input: all cubes need to have the " +
                                 "same dimensions")
        self.y = np.ravel(data1.data)

        # Matplotlib
        plt.style.use(MPLSTYLE)
        self.fig, self.ax = plt.subplots()

###############################################################################

    def plot(self, x_label=None, y_label=None, title=None, add_info=None):
        """
        Arguments
            x_label  : label of x-axis
            y_label  : label of y-axis
            title    : title of the plot
            add_info : dictionay with additional text

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180209-A_schl_ma: written
        """

        # Parse arguments
        if (x_label is None):
            x_label = 'data 2'
        if (y_label is None):
            y_label = 'data 1'

        # Create scatter plot
        for d in self.data2:
            x = np.ravel(d.data)
            self.ax.plot(x, self.y, linestyle='none',
                         marker='o', markerfacecolor='none')

        # Plot appearance
        self.ax.set_xlabel(x_label)
        self.ax.set_ylabel(y_label)
        if (title is not None):
            self.ax.set_title(title)
        if (add_info is not None):
            text = ""
            for key in add_info:
                text += "{0} = {1}\n".format(key, add_info[key])
            self.ax.text(0.02, 0.8, text, transform=self.ax.transAxes)

        self.fig.tight_layout()
        return self.fig


class Plot2D_deprecated(object):
    """
    Description
        Basic class for 2-dimensional plotting

    Contents
        method plot
    """

    # Class attributes
    LATS = ['latitude']     # accepted lat names
    LONS = ['longitude']    # accepted lon names
    TIME = ['time']         # accepted time names
    GRIDSPEC_LENGTH = 5     # size of the matplotlib gridspec

###############################################################################

    def __init__(self, cube):
        """
        Arguments
            cube : iris cube

        Description
            Initializes the class.

        TODO
            optional summary placement
            swap x/y axis

        Modification history
            20180207-A_schl_ma: written
        """

        # Check arguments
        if (not isinstance(cube, iris.cube.Cube)):
            raise TypeError("Invalid input: expected iris cube")
        self.cube = iris.util.squeeze(cube)
        if (self.cube.ndim != 2):
            raise TypeError("Invalid input: expected 2-dimensional iris cube")
        try:
            self.name = cube.long_name
        except:
            pass
        self.units = ' [' + str(cube.units) + ']'


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
            self.summary_location = 'right'

        # Lat/time plot
        elif (self.lat_var is not None and self.time_var is not None):
            self.plot_type = 'lattime'
            self.summary_location = 'right'

        # Lon/time plot
        elif (self.lon_var is not None and self.time_var is not None):
            self.plot_type = 'lontime'
            self.summary_location = 'bottom'

        # Default case
        else:
            raise TypeError("Invalid input: cube does not contain supported " +
                            "dimensions")

        # Setup matplotlib
        plt.style.use(MPLSTYLE)

###############################################################################

    def plot(self, summary_plot=False, colorbar_ticks=None, x_label=None,
             y_label=None, title=None, ax=None, fig=None):
        """
        Arguments
            summary_plot   : Add summary line plot
            colorbar_ticks : Ticks of the colorbar
            x_label        : label of x-axis
            y_label        : label of y-axis
            title          : title of the plot

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180207-A_schl_ma: written
        """


        if ax is None:

            # Summary plot yes/no
            if (not summary_plot):
                G = gridspec.GridSpec(1, 1)
                colorbar = 'vertical'
            elif (self.summary_location == 'right'):
                colorbar = 'horizontal'
                G = gridspec.GridSpec(2, 2, width_ratios=[4, 1],
                                      height_ratios=[8, 1])
            elif (self.summary_location == 'bottom'):
                colorbar = 'vertical'
                G = gridspec.GridSpec(2, 2, width_ratios=[8, 1],
                                      height_ratios=[4, 1])

            # Main plot
            plt.subplot(G[0])
            if (title is not None):
                plt.title(title)
            else:
                plt.title(self.name + self.units)
            if (x_label is not None):
                plt.xlabel(x_label)
            if (y_label is not None):
                plt.ylabel(y_label)
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
            if (summary_plot):
                ax_cb = plt.subplot(G[cb_axes])
                cb = plt.colorbar(ax_main, cax=ax_cb, orientation=colorbar,
                                  ticks=colorbar_ticks)
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
                cb = plt.colorbar(ax_main, orientation=colorbar,
                                  ticks=colorbar_ticks)

            plt.tight_layout()
            return plt.gcf()

        else:

            if fig is None:
                assert False, "fig should be provided and is missing!"

            # Summary plot yes/no
            if (not summary_plot):
                G = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=ax)
                colorbar = 'vertical'
            elif (self.summary_location == 'right'):
                colorbar = 'horizontal'
                G = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=ax, width_ratios=[4, 1],
                                      height_ratios=[8, 1])
            elif (self.summary_location == 'bottom'):
                colorbar = 'vertical'
                G = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=ax, width_ratios=[8, 1],
                                      height_ratios=[4, 1])
            # Main plot
            axl = plt.Subplot(fig, G[0])
            if (title is not None):
                plt.title(title)
            else:
                plt.title(self.name + self.units)
            if (x_label is not None):
                plt.xlabel(x_label)
            if (y_label is not None):
                plt.ylabel(y_label)
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
            fig.add_subplot(ax_main)

            # Line plot
            if (summary_plot):
                ax_cb = plt.Subplot(fig, G[cb_axes])
                cb = plt.colorbar(ax_main, cax=ax_cb, orientation=colorbar,
                                  ticks=colorbar_ticks)
                fig.add_subplot(ax_cb)
                axl = plt.Subplot(fig, G[line_plot_axes])
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
                fig.add_subplot(axl)
            else:
                cb = plt.colorbar(ax_main, orientation=colorbar,
                                  ticks=colorbar_ticks)
                fig.add_subplot(cb)

            return


class Plot2D(object):
    """
    Description
        Basic class for 2-dimensional plotting

    Contents
        method plot
    """

    # Class attributes
    LATS = ['latitude']     # accepted lat names
    LONS = ['longitude']    # accepted lon names
    TIME = ['time']         # accepted time names
    LEVS = ['air_pressure'] # accepted level names
    MAX_COLUMNS = 3
    TIME_LABEL = 'Date'
    TIME_FORMAT = '%Y-%m-%d'


    def __init__(self, cubes):
        """
        Arguments
            cubes : iris cube or a list of iris cubes

        Description
            Initializes the class.

        TODO
            optional summary placement
            swap x/y axis

        Modification history
            20180207-A_muel_bn: copied Plot2D and adjusted
        """

        # Check arguments
        try:
            self.n_cubes = len(cubes)
        except TypeError:
            self.n_cubes = 1
            cubes = [cubes]
        for (idx, cube) in enumerate(cubes):
            if not isinstance(cube, iris.cube.Cube):
                raise TypeError("Invalid input: expected single iris cube or "
                                "list of iris cubes")
            cubes[idx] = iris.util.squeeze(cube)
            if cubes[idx].ndim != 2:
                raise TypeError("Invalid input: expected 2-dimensional iris "
                                "cubes at list index {}".format(idx))
        self.cubes = cubes
        if self.n_cubes < 1:
            raise TypeError("Invalid input: expected at least one cube, not "
                            "an empty list")

        # Get properties of the cubes
        self.names = [None] * self.n_cubes
        self.units = [None] * self.n_cubes
        self.time_vars = [None] * self.n_cubes
        self.lev_vars = [None] * self.n_cubes
        self.lat_vars = [None] * self.n_cubes
        self.lon_vars = [None] * self.n_cubes
        for (idx, cube) in enumerate(self.cubes):

            # Name
            try:
                self.names[idx] = cube.long_name
            except:
                pass

            # Units
            self.units[idx] = (' [' + str(cube.units) + ']')

            # Dimension names
            dim_names = [dim.standard_name for dim in cube.dim_coords]
            for dim in dim_names:
                if (dim in self.__class__.TIME):
                    self.time_vars[idx] = dim
                elif (dim in self.__class__.LEVS):
                    self.lev_vars[idx] = dim
                elif (dim in self.__class__.LATS):
                    self.lat_vars[idx] = dim
                elif (dim in self.__class__.LONS):
                    self.lon_vars[idx] = dim

        # Lat/lon plot
        if (all(lat is not None for lat in self.lat_vars) and
                all(lon is not None for lon in self.lon_vars)):
            self.plot_type = 'latlon'

        # Lat/time plot
        elif (all(lat is not None for lat in self.lat_vars) and
              all(time is not None for time in self.time_vars)):
            self.plot_type = 'lattime'

        # Lon/time plot
        elif (all(lon is not None for lon in self.lon_vars) and
              all(time is not None for time in self.time_vars)):
            self.plot_type = 'lontime'

        # Lat/lev plot
        elif (all(lat is not None for lat in self.lat_vars) and
              all(lev is not None for lev in self.lev_vars)):
            self.plot_type = 'latlev'

        # Lon/lev plot
        elif (all(lon is not None for lon in self.lon_vars) and
              all(lev is not None for lev in self.lev_vars)):
            self.plot_type = 'lonlev'

        # time/lev plot
        elif (all(time is not None for time in self.time_vars) and
              all(lev is not None for lev in self.lev_vars)):
            self.plot_type = 'timelev'

        # Default case
        else:
            raise TypeError("Invalid input: Not all cubes have the same "
                            "dimensions or at least on of them has an "
                            "unknown dimension name")

        # Setup matplotlib
        plt.style.use(MPLSTYLE)

###############################################################################

    def plot(self, summary_plot=False, colorbar_ticks=None, x_label=None,
             y_label=None, title=None, ax=None, fig=None, vminmax=None,
             color=None, color_type=None, color_reverse=False):
        """
        Arguments
            summary_plot   : Add summary line plot
            colorbar_ticks : Ticks of the colorbar
            x_label        : label of x-axis
            y_label        : label of y-axis
            title          : title of the plot

        Returns
            Matplotlib figure

        Description
            Actual plotting routine

        Modification history
            20180515-A_muel_bn: copied Plot2D and adjusted
        """

        # Preprocessing cube information
        if self.n_cubes == 1:
            self.cubes[0].rename(title)
        for (idx, cube) in enumerate(self.cubes):
            try:
                self.cubes[idx] = cube.intersection(longitude=(-180, 180))
            except:
                pass

        # Color
        if color is None:
            brewer_cmap = mpl_cm.get_cmap('brewer_Spectral_11')
        else:
            if color_type is None or color_type not in color.keys():
                if color_reverse:
                    col_save=color["default"]
                    color["default"]=color["default"]+"_r"
                brewer_cmap = mpl_cm.get_cmap(color["default"], lut=11)
            else:
                if color_reverse:
                    col_save=color[color_type]
                    color[color_type]=color[color_type]+"_r"
                brewer_cmap = mpl_cm.get_cmap(color[color_type], lut=11)

        # Vminmax
        if vminmax is None:
            vmin, vmax = 1E99, -1E99
            for cube in self.cubes:
                try:
                    temp_min, temp_max = np.nanpercentile(cube.data.data[
                        np.logical_not(cube.data.mask)], [5.0, 95.0])
                    if (temp_min < vmin):
                        vmin = temp_min
                    if (temp_max > vmax):
                        vmax = temp_max
                except:
                    pass
            if vmin > vmax:
                vmin = None
                vmax = None
                levels = None
        else:
            vmin, vmax = vminmax
        if vmin is not None:
            rounder = int(np.ceil(-np.log10(vmax - vmin) + 1))
            vmin, vmax = np.round([vmin, vmax], rounder)
            levels = np.round(np.linspace(vmin, vmax, num=11), rounder)

        # Check axes
        if ax is None:
            ax = [plt.gca()]
        try:
            if len(ax) == 2:
                summary_plot = True
            elif len(ax) > 2:
                raise ValueError("Invalid input: axes should not be more "
                                 "than 2!")
        except TypeError:
            ax = [ax]

        # Plot summary if necessary
        if summary_plot:
            if len(ax) < 2:
                raise ValueError("Invalid input: Need 2 axes for summary plot")
            plt.sca(ax[1])

            # Latlon plot and multiple cubes not supported yet
            if self.n_cubes > 1:
                raise ValueError("Invalid input: summary plot for multiple "
                                 "cubes is not supported")
            cube = self.cubes[0]
            time_var = self.time_vars[0]
            lev_var = self.lev_vars[0]
            lat_var = self.lat_vars[0]
            lon_var = self.lon_vars[0]

            # Get variable to collapse
            if (self.plot_type == 'latlon'):
                raise ValueError("Invalid input: latlon should not have a "
                                 "summary plot")
            elif (self.plot_type == 'lattime'):
                collapse = time_var
            elif (self.plot_type == 'lontime'):
                collapse = time_var
            else:
                collapse = [c.name() for c in cube.coords() if c.name() !=
                            lev_var and len(c.points) > 1][0]
            if (collapse == lat_var):
                grid_areas = iris.analysis.cartography.area_weights(cube)
            else:
                grid_areas = None
            cube_line = cube.collapsed(collapse, iris.analysis.MEAN,
                                       weights=grid_areas)

            # Lattime
            if self.plot_type == 'lattime':
                x = cube_line
                y = cube_line.coord(lat_var)
                latrange = cube.coords(lat_var).pop()
                # lat_inc = np.diff(latrange.points).mean()
                latrange = (np.min(latrange.points), np.max(latrange.points))
                plt.gca().set_xlim(vmin, vmax)
                plt.gca().set_ylim(latrange)

            # Lontime
            elif self.plot_type == 'lontime':
                x = cube_line.coord(lon_var)
                y = cube_line
                lonrange = cube.coords(lon_var).pop()
                # lon_inc = np.diff(lonrange.points).mean()
                lonrange = (np.min(lonrange.points), np.max(lonrange.points))
                plt.gca().set_ylim(vmin, vmax)
                plt.gca().set_xlim(lonrange)
            # Lev plots
            else:
                x = cube_line
                y = cube_line.coord(lev_var)
                levrange = cube.coords(lev_var).pop()
                # lat_inc = np.diff(latrange.points).mean()
                levrange = (np.min(levrange.points), np.max(levrange.points))
                plt.gca().set_xlim(vmin, vmax)
                plt.gca().set_ylim(levrange)

            # Plot
            iplt.plot(x, y)

        # Setup axes for plotting multiple cubes
        n_columns = min([self.n_cubes, self.__class__.MAX_COLUMNS])
        n_rows = (self.n_cubes - 1) // self.__class__.MAX_COLUMNS + 1
        gs = gridspec.GridSpec(n_rows, n_columns)

        # Iterate over cubes
        for (idx, cube) in enumerate(self.cubes):
            column = idx % self.__class__.MAX_COLUMNS
            row = idx // self.__class__.MAX_COLUMNS
            if self.n_cubes > 1:
                plt.sca(plt.subplot(gs[row, column]))
            else:
                plt.sca(ax[0])

            # Plot map
            # (this needs to be done due to an error in cartopy)
            try:
                # qplt.pcolormesh(self.cubes, cmap=brewer_cmap, vmin=vmin,
                #                 vmax=vmax) #, levels=levels, extend='both')
                # qplt.contourf(self.cubes, cmap=brewer_cmap, vmin=vmin,
                #               vmax=vmax, levels=levels, extend='both')
                qplt.pcolormesh(cube, cmap=brewer_cmap, vmin=vmin,
                                vmax=vmax)
                if 'time' in self.plot_type:
                    time_coord = cube.coord(self.time_vars[idx])
                    if self.plot_type == 'lontime':
                        (locs, _) = plt.yticks()
                    else:
                        (locs, _) = plt.xticks()
                    labels = unit.num2date(locs,
                                           time_coord.units.name,
                                           time_coord.units.calendar)
                    print labels
#                    for (idx, _) in enumerate(labels):
#                        labels[idx] = labels[idx].strftime(
#                            self.__class__.TIME_FORMAT)
#                    if self.plot_type == 'lontime':
#                        plt.yticks(locs, labels)
#                        plt.ylabel(self.__class__.TIME_LABEL)
#                    else:
#                        plt.xticks(locs, labels)
#                        plt.xlabel(self.__class__.TIME_LABEL)
            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                print(exc_type, fname, exc_tb.tb_lineno)
                qplt.pcolormesh(cube, cmap=brewer_cmap, vmin=vmin, vmax=vmax)
                plt.text(0.5, 0.5,'Data cannot be displayed as intended due '
                                  'to cartopy bug! \n Deviations are color '
                                  'levels and time axis display. \n Future '
                                  'updates of cartopy module may resolve '
                                  'this issue (#946).',
                                  horizontalalignment='center',
                                  verticalalignment='center',
                                  transform=plt.gca().transAxes)
            if self.plot_type == 'latlon':
                mean = cube.collapsed([coord.name() for coord in
                                       cube.coords()],
                                      iris.analysis.MEAN).data
                std = cube.collapsed([coord.name() for coord in cube.coords()],
                                     iris.analysis.STD_DEV).data
                plt.gca().coastlines()
                plt.gca().gridlines(crs=ccrs.Geodetic(), color="k",
                                    linestyle=':')
                plt.gca().text(-180, -100,
                               r'mean: {0} $\pm$ {1} '.format(
                                   str(round(mean, 2)),str(round(std, 2))))

            # Label plots (supports at most 26 figures at the moment)
            if self.n_cubes > 1:
                letter = '(' + ascii_lowercase[idx % 26] + ')'
                plt.text(0.01, 0.99, letter,
                        horizontalalignment='left',
                        verticalalignment='top',
                        transform=plt.gca().transAxes)

        # Colors
        if color_type is None or color_type not in color.keys():
            if color_reverse:
                color["default"] = col_save
        else:
            if color_reverse:
                color[color_type] = col_save

        # Set position of summary plot
        if summary_plot:
            if self.plot_type not in ['lontime', 'latlon']:
                bb = ax[1].get_position()
                bb.y0 = ax[0].get_position().y0
                ax[1].set_position(bb)

        # Tight layout
        try:
            plt.tight_layout()
        except: 
            pass
        return


class Plot2D_blank(Plot2D):
    """
    Description
        Blank class for 2-dimensional plotting

    Contents
        method plot
    """
    def __init__(self, cubes):
        """
        Arguments
            cube : iris cube

        Description
            Initializes the class.

        TODO
            optional summary placement
            swap x/y axis

        Modification history
            20180515-A_muel_bn: written
        """
        super(Plot2D_blank, self).__init__(cubes)
        # erase all data
        for cube in self.cubes:
            cube.data=cube.data*np.nan


class Plot1D(object):
    """
    Description
        Basic class for 21-dimensional plotting

    Contents
        method plot
    """

    # Class attributes
    LATS = ['latitude']     # accepted lat names
    LONS = ['longitude']    # accepted lon names
    TIME = ['time']         # accepted time names


    def __init__(self, cube):
        """
        Arguments
            cube : iris cube

        Description
            Initializes the class.


        Modification history
            20180527-A_muel_bn: copied Plot2D_2 and adjusted
        """

        # Check arguments
        if (not isinstance(cube, iris.cube.Cube)):
            raise TypeError("Invalid input: expected iris cube")
        self.cube = iris.util.squeeze(cube)
        if (self.cube.ndim != 1):
            raise TypeError("Invalid input: expected 1-dimensional iris cube")
        try:
            self.name = cube.long_name
        except:
            pass
        self.units = ' [' + str(cube.units) + ']'

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
        if (self.lat_var is not None):
            self.plot_type = 'lat'

        # Lat/time plot
        elif (self.time_var is not None):
            self.plot_type = 'time'

        # Lon/time plot
        elif (self.lon_var is not None):
            self.plot_type = 'lon'

        # Default case
        else:
            raise TypeError("Invalid input: cube does not contain supported " +
                            "dimensions")

        # Setup matplotlib
        plt.style.use(MPLSTYLE)

###############################################################################

    def plot(self, summary_plot=False, colorbar_ticks=None, x_label=None,
             y_label=None, title=None, ax=None, fig=None, vminmax=None):
        """
        Arguments
            summary_plot   : Add summary line plot
            colorbar_ticks : Ticks of the colorbar
            x_label        : label of x-axis
            y_label        : label of y-axis
            title          : title of the plot

        Returns
            Matplotlib figure instance

        Description
            Actual plotting routine

        Modification history
            20180527-A_muel_bn: copied Plot2D_2 and adjusted
        """

        # preprocessing cube information
        self.cube.rename(title)

        #brewer_cmap = mpl_cm.get_cmap('brewer_Spectral_11')

        # check axes
        if ax is None:
            ax = [plt.gca()]
        try:
            if len(ax) >= 2:
                raise ValueError("Invalid input: axes should not be more "
                                 "than 1!")
        except TypeError:
            ax = [ax]

        plt.sca(ax[0])

        # plot line
        try:
#            print self.cube
#            print self.cube.data
#            print self.cube.coords("time")[0].points
            plt.plot(self.cube.coords("time")[0].points,self.cube.data)
            plt.title(title)
            plt.grid()

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print(exc_type, fname, exc_tb.tb_lineno)
            print 'We did not expect this to fail!'
            plt.plot(self.cube.coords("time").points,self.cube.data)
        plt.tight_layout()

        return
