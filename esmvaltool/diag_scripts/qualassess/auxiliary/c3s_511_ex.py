#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:59:59 2018

@author: bmueller
"""

import iris
import iris.pandas as ipd
import os
import sys
import matplotlib.pyplot as plt
import datetime
import numpy as np
import cf_units
from scipy import sin, cos, tan, arctan, arctan2, arccos, pi, deg2rad

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .libs.c3s_511_util import cube_sorted
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup

class ex_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(ex_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
        self.__regions__ = dict({
            'CE_JJA_2003': {  # https://en.wikipedia.org/wiki/2003_European_heat_wave
                'latitude': (43, 47),
                'longitude': (-2, 4),
                'time': (datetime.datetime(2003, 7, 20),
                         datetime.datetime(2003, 8, 20)
                         )
                }})
    
    def run_diagnostic(self):
#        self.sp_data = self.__spatiotemp_subsets__(self.sp_data)['Europe_2000']
        self.__do_extremes__()

        self.__do_full_report__()

    def __do_extremes__(self):
        
        this_function =  "extremes example"
        
        which_percentile = 10
        window_size = 5 # one directional 5 => 11
        masked_val = None
        
        # this the extremes example
        
        # handling of the cube restrictions copied from the basic diagnostic
        # needs to be set locally
        #################
        if not self.var3D:
            cube = self.sp_data
        else:
            cube = self.sp_data.extract(iris.Constraint(
                    coord_values={str(self.level_dim): lambda cell: cell == max(self.levels)}))

        # adjustment to ids and filenames
        if self.var3D:
            basic_filename = self.__basic_filename__ + "_lev" + str(max(self.levels))
            dataset_id = [self.__dataset_id__[0], "at", "level", str(
                max(self.levels)), str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        ##################
        
        # we assume that the regions fully cover the event in space and time (TODO: build in a check for this)
        list_of_plots = []
        
        # Loop over the different regions (i.e. the different events)
        for r,def_r in self.__regions__.items():
            # Now define the three cubes. Note that now they are really cubes,
            # not dictionaries that contain cubes.
            event_cube = self.__spatiotemp_subsets__(cube,{r:def_r})[r]

            spatial_def = def_r.copy()
            spatial_def.pop("time")
            clim_cube = self.__spatiotemp_subsets__(cube,{r:spatial_def})[r]
            ex_cube = event_cube.copy()

            # Select this specific region, note that loc_cube is a dictionary with loc_cube[key]
            # holding the actual cubes.
            loc_cube = self.__spatiotemp_subsets__(cube,{r:def_r})

            
            list_of_cubes = []
            
            # Now loop over each gridpoint in the selected region
            counter_gridpoints = 0
            n_gridpoints = event_cube.shape[1]*event_cube.shape[2]
            self.__logger__.info("Start calculation of extreme climatology \
                                 for %s gridpoints",n_gridpoints)

            for ii in range(event_cube.shape[1]):
                for jj in range(event_cube.shape[2]):
                    # Get doy
                    gridpoint_ts = ipd.as_series(event_cube[:,ii,jj])
                    doy_to_process = event_cube.coord('day_of_year').points
                    n_doy = len(doy_to_process)
                    for n,doy in enumerate(event_cube.coord('day_of_year').points):
                        doy_window = (doy - window_size < gridpoint_ts.index.dayofyear) &\
                                     (gridpoint_ts.index.dayofyear < doy + window_size)
                        perc_val = np.nanpercentile(gridpoint_ts[doy_window],which_percentile)
                        ex_cube.data[n,ii,jj] = perc_val
                    self.__logger__.info("Progress of xclim: %s percent",\
                                         np.round(100.*(counter_gridpoints/n_gridpoints),decimals=1))
                    counter_gridpoints += 1
            #TODO think about propagation of nan values
            ex_cube.data = np.ma.masked_equal(ex_cube.core_data(), np.nan)
            self.__logger__.info("Finished calculation of extreme climatology")

            # Calculate exceedance of the extremes climatology
            if which_percentile > 50: # we are interested in over threshold values
                incident_cube = event_cube - ex_cube
            elif which_percentile < 50: # we are interested in under threshold values
                incident_cube = ex_cube - event_cube
            else:
                self.__loger__.error("Percentile can not be 50, that wouldn't be an extreme climatology")
                raise ValueError
            # Note that due to the above check, incident_cube values of the event are always positive
            # therefore mask negative values
            incident_cube.data = np.ma.masked_where(incident_cube.core_data() <= 0, incident_cube.core_data(), copy = True)
            
            # Sidenote:
            # TODO: Insert a check that incident_cube does not have two or more dims with the same shape
            # otherwise broadcast could fail and produce erroneous results without notice
            
            # calculate severity
            severity = incident_cube * np.atleast_3d(np.array([np.diff(bds) for bds in loc_cube[r].coord("time").bounds]))
            severity.units = cf_units.Unit(str(cube.units) + 
                                           " " + 
                                           str(cube.coord("time").units).split(" ")[0])
            severity = severity.collapsed("day_of_year", iris.analysis.SUM) 
            severity.long_name = "severity"
            
            # calculate magnitude
            magnitude = incident_cube.collapsed("day_of_year", iris.analysis.MAX)
            magnitude.long_name = "magnitude"
            
            # calculate duration
            duration = (incident_cube * 0 + 1) * np.atleast_3d(np.array([np.diff(bds) for bds in loc_cube[r].coord("time").bounds]))
            duration.units = (str(cube.coord("time").units).split(" ")[0])
            duration = duration.collapsed("day_of_year", iris.analysis.SUM)
            duration.long_name = "duration"
            
            # calculate spatial averages
            grid_areas = iris.analysis.cartography.area_weights(severity)
            #TODO: below three lines should be changed to 
            # a nanmean (i.e. taking the mean over the not-nan values)
            severity_av = severity.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas)             
            magnitude_av = magnitude.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas) 
            duration_av = duration.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas) 
            
            # calculate extent
            extent = ((duration * 0 + 1.) * grid_areas).collapsed(["latitude","longitude"], iris.analysis.SUM)/1e6
            extent.units = cf_units.Unit("km2")
            
            # set up table
            # TODO export csv
            self.__logger__.info("Extremes table")
            self.__logger__.info("mean severity: {:.2f} {}".format(severity_av.data, severity_av.units))
            self.__logger__.info("mean magnitude: {:.2f} {}".format(magnitude_av.data, magnitude_av.units))
            self.__logger__.info("mean duration: {:.2f} {}".format(duration_av.data, duration_av.units))
            self.__logger__.info("extent: {:.2f} {}".format(extent.data, extent.units))
            # plotting for trials
            import iris.quickplot as qplt
            
            for dat in ["severity", "magnitude", "duration"]:
                #self.__logger__.info(locals()[dat])
            
                qplt.pcolormesh(locals()[dat])#, vmin = 250, vmax = 300)
                
                # Add coastlines to the map created by contourf.
                plt.gca().coastlines("10m")
                plt.colorbar()
                plt.savefig(self.__plot_dir__ + os.sep + dat + ".png")
                
                plt.close()

            import IPython;IPython.embed()

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_extremes_input",
                                      expfile="_extremes.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:

            self.reporting_structure.update(
                    {"Extremes": 
                        {"plots": list_of_plots,
                         "freetext": expected_input}})
        else:

            self.reporting_structure.update(
                    {"Extremes": 
                        {"plots": list_of_plots}})
