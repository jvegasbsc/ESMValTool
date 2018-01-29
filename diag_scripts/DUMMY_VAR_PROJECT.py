# -*- coding: utf-8 -*-

"""
###############################################################################
DUMMY_VAR_PROJECT.py
Author: Birgit Hassler (DLR, Germany)
C3S_511 service
###############################################################################

Description
    Test file for basic diagnostics for a single ESC.

Required diag_script_info attributes (diagnostics specific)
    [ecs_plots]
        plot : Switch to plot the linear regression needed for the ECS
               calculation
    [netcdf]
        filename  : Name of the output file
        overwrite : Overwrite existing files

Optional diag_script_info attributes (diagnostic specific)
    [main_plot]
        fontsize : Fonzsize used in the plot
        xmin     : Left boundary of the plot
        xmax     : Right boundary of the plot
        ymin     : Lower boundary of the plot
        ymax     : Upper boundary of the plot
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
    20180125-A_hass_bg: written

###############################################################################
"""


# ESMValTool python packages
from auxiliary import info, warning, error
from esmval_lib import ESMValProject
from ESMValMD import ESMValMD
from grid_operations import GridOperations

# NetCDF4
from netCDF4 import Dataset

# Basic python packages
from collections import OrderedDict
from datetime import datetime
from scipy import stats
import ConfigParser
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable




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

	TOZ = "toz"

	VARIABLES = [(TOZ)]


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
		warning("Selected file type for plots is not supported",
                verbosity, 0, exit_on_warning)
		plot_file_type = "ps"

    # Check if all needed varibles are present
	for var in VARIABLES:
		if (var not in vars):
			error("no data for variable '{0}' available, ".format(var) + \
                  "please check your namelist")

    # Write references:
	E.write_references(diag_name,                # diagnostic script name
						["A_hass_bg"],           # authors
						[""],                    # contributors
						["D_0002"],              # diagnostic
						[""],                    # observations
						["P_square4ECVs"],       # project
						project_info,
						verbosity,
						False)

	# Read configuration file
	modelconfig = ConfigParser.ConfigParser()
	modelconfig.read(config_file)

    # Get all models (in our case the observations)
	models = E.get_all_clim_models([var for var in VARIABLES])
	print(models)
	
	###########################################################################
    # Collect data of the models
    ###########################################################################

	info("", verbosity, 1)
	info("Starting calculation", verbosity, 1)
	info("", verbosity, 1)

	#tco_data = OrderedDict((exp, OrderedDict()) for exp in VARIABLES[TOZ])
	units = OrderedDict()
	climo_files = []

    # Try to get necessary model information
	for model_path in models:
		E.add_to_filelist(model_path)
		climo_files.append(model_path)
		model_info = models[model_path]

        # Try to get necessary model information
		model_var = None
		model_name = None
		model_exp = None
		try:
			model_var = model_info["var"]
			model_name = model_info["name"]
		except KeyError:
			warning("Could not retrieve all desired model information of " +
                    "model {0}".format(model_info), verbosity, 0,
                    exit_on_warning)
			continue

		# Skip unnecesseary calculations
		grid_op = GridOperations(model_path, model_var, global_conf)
		info("Retrieving '{0}' from model '{1}' [{2}]".format(model_var,
                                                         model_name,
														 model_exp),
			 verbosity, 1)
		units.update({model_var: grid_op.get_var_units()})
		
		# Get netCDF latitudes and longitudes
		#print(model_path)
		netCDF_data = Dataset(model_path, 'r')
		#print netCDF_data.variables['lat']
		
		lats = netCDF_data.variables['lat'][:]
		lons = netCDF_data.variables['lon'][:]
		data = netCDF_data.variables['toz'][:]
		#lons_idx = np.where(lons > 180)
		lons[lons > 180] =  lons[lons > 180] - 360.
		lons_adj = np.concatenate((lons[180:],lons[:180]), axis = 0)
		#print(lats)
		#print(lons_adj)
		#print(data.shape)

        # mean
		if (model_var == TOZ):
			tco_mean = grid_op.average(spatial_axis="all", period="annual",
									   spatial_weighting=True, region="global")

			#tco_grid_mean = grid_op.average(period="total",
			#						        spatial_weighting=True, region="global")
			tco_grid_mean = np.mean(data, axis=0)
			#print(tco_grid_mean.shape)
			tco_grid_mean_adj = np.concatenate((tco_grid_mean[:,180:],tco_grid_mean[:,:180]), axis = 1)
			#print(tco_grid_mean_adj.shape)
			tco_grid_std = np.std(data, axis=0)
			tco_grid_std_adj = np.concatenate((tco_grid_std[:,180:],tco_grid_std[:,:180]), axis = 1)

			tco_grid_num = sum(~np.isnan(data))
			tco_grid_num_adj = np.concatenate((tco_grid_num[:,180:],tco_grid_num[:,:180]), axis = 1)
			print(tco_grid_num.shape)
			
    # Empty line
	info("", verbosity, 1)

	# Matplotlib instance
	fig, axes = plt.subplots()

	# Plot line data
	#----------------------
	if (write_plots):
		info("Create global mean TCO plot '{0}'".format(plot_dir),
             verbosity, 1)
	
        # Get values from configuration file
		main_plot_section = "main_plot"
		cfg_options = {"filename": "tco_mean", "fontsize": 18.0,
                       "xmin": 260.0, "xmax": 280.0,
                       "ymin": 1995.0, "ymax": 2012.0}
		cfg = E.get_config_options(modelconfig, main_plot_section, cfg_options)

		# Line plot
		years = []
		for i in range(2010 - 1997 + 1):
			years.append(1997 + i)
		
		#print(years)
		
		axes.plot(years, tco_mean,
                  linestyle = "solid", linewidth = 5, color = "red")
		
		
        # Save line plot
		filename = cfg["filename"] + "." + plot_file_type
		fig.savefig(plot_dir + filename, 
                    bbox_inches="tight", orientation="landscape")
		plt.close("all")
		
		# Plot map (mean)
		#-----------------------------------
		#-- create figure and axes instances
		dpi = 100
		fig = plt.figure(figsize=(1100/dpi, 1100/dpi), dpi=dpi)
		ax  = fig.add_axes([0.1,0.1,0.8,0.9])

		#-- create map
		map = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
              resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
		
		#-- draw coastlines, state and country boundaries, edge of map
		map.drawcoastlines()
		#map.drawstates()
		#map.drawcountries()

		#-- create and draw meridians and parallels grid lines
		map.drawparallels(np.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=10)
		map.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10)
		
		#-- plot the contours on the map
		#print(lons.shape)
		#print(lats.shape)
		#print(np.transpose(tco_grid_mean).shape)
		cnplot = plt.contourf(lons_adj, lats, tco_grid_mean_adj, cmap=cm.jet)

		#-- add plot title
		plt.title('Mean TCO')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="3%", pad=0.05)
		# Make a colorbar for the ContourSet returned by the contourf call.
		cbar = plt.colorbar(cnplot, cax=cax)
		cbar.ax.set_ylabel('TCO [DU]')

		# Save map
		filename = "TCO_mean_map" + "." + plot_file_type
		plt.savefig(plot_dir + filename, bbox_inches="tight", dpi=dpi)
		plt.close("all")

		
		# Plot map (std)
		#-----------------------------------
		#-- create figure and axes instances
		dpi = 100
		fig = plt.figure(figsize=(1100/dpi, 1100/dpi), dpi=dpi)
		ax  = fig.add_axes([0.1,0.1,0.8,0.9])

		#-- create map
		map = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
              resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
		
		#-- draw coastlines, state and country boundaries, edge of map
		map.drawcoastlines()
		#map.drawstates()
		#map.drawcountries()

		#-- create and draw meridians and parallels grid lines
		map.drawparallels(np.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=10)
		map.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10)
		
		#-- plot the contours on the map
		cnplot = plt.contourf(lons_adj, lats, tco_grid_std_adj, cmap=cm.jet)

		#-- add plot title
		plt.title('Standard deviation of TCO')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="3%", pad=0.05)
		# Make a colorbar for the ContourSet returned by the contourf call.
		cbar = plt.colorbar(cnplot, cax=cax)
		cbar.ax.set_ylabel('TCO [DU]')

		
		# Save map
		filename = "TCO_std_map" + "." + plot_file_type
		plt.savefig(plot_dir + filename, bbox_inches="tight", dpi=dpi)
		plt.close("all")


		# Plot map (num)
		#-----------------------------------
		#-- create figure and axes instances
		dpi = 100
		fig = plt.figure(figsize=(1100/dpi, 1100/dpi), dpi=dpi)
		ax  = fig.add_axes([0.1,0.1,0.8,0.9])

		#-- create map
		map = Basemap(projection='cyl',llcrnrlat= -90.,urcrnrlat= 90.,\
              resolution='c',  llcrnrlon=-180.,urcrnrlon=180.)
		
		#-- draw coastlines, state and country boundaries, edge of map
		map.drawcoastlines()
		#map.drawstates()
		#map.drawcountries()

		#-- create and draw meridians and parallels grid lines
		map.drawparallels(np.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=10)
		map.drawmeridians(np.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=10)
		
		#-- plot the contours on the map
		cnplot = plt.contourf(lons_adj, lats, tco_grid_num_adj, cmap=cm.jet)

		#-- add plot title
		plt.title('Number of available data points')

		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="3%", pad=0.05)
		# Make a colorbar for the ContourSet returned by the contourf call.
		cbar = plt.colorbar(cnplot, cax=cax)
		cbar.ax.set_ylabel('number')

		
		# Save map
		filename = "TCO_num_map" + "." + plot_file_type
		plt.savefig(plot_dir + filename, bbox_inches="tight", dpi=dpi)
		plt.close("all")
		
	print("test successful")