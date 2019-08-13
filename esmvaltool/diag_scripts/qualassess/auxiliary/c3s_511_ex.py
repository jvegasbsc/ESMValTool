#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:59:59 2018

@authors: bmueller, bcrezee
"""

import pandas as pd
import iris
import iris.pandas as ipd
import iris.quickplot as qplt
import os
import re
import sys
import matplotlib.pyplot as plt
import datetime
import numpy as np
import cf_units
from scipy import sin, cos, tan, arctan, arctan2, arccos, pi, deg2rad

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .libs.c3s_511_util import cube_sorted, read_extreme_event_catalogue
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup
    
from multiprocessing import Pool
import itertools as it


class ex_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(ex_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
        
        # all required input can be extracted from the extremes dictionary
        self.__logger__.info(self.__extremes__)
        
        # Initialize extremes regions as empty
        self.__extremes_regions__ = dict()

    def run_diagnostic(self):
#        self.sp_data = self.__spatiotemp_subsets__(self.sp_data)['Europe_2000']
        self.__do_extremes__()

        self.__do_full_report__()

    def __do_extremes__(self):
        
        this_function =  "extremes example"

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
        
        list_of_plots = []
        
        # end of standard handling
                
        # Read settings from recipe
        min_measurements = self.__extremes__["min_measurements"]
        which_percentile = self.__extremes__["which_percentile"]
        window_size = self.__extremes__["window_size"]
        extreme_events = self.__extremes__["extreme_events"]
        num_processors = self.__extremes__["multiprocessing"]
        
        self.__logger__.info("Reading extreme event table")
        ex_table,raw_table = read_extreme_event_catalogue()
        self.__logger__.info("Finished parsing extreme event table")

        # Loop through the events
        self.__logger__.info("Adding selected events to list for processing: ")
        # Now start adding the events to the dictionary for further processing
        if type(extreme_events) is not list:
            extreme_events = [extreme_events]
        for extreme_event_id in extreme_events:
            self.__logger__.info("%s",extreme_event_id)
            try:
                single_event = ex_table.loc[extreme_event_id].to_dict()
                # Now prepare dictionary to be added to regions
                single_event_as_region = {}
                single_event_as_region['latitude'] = (single_event['Lat_from'],\
                                                      single_event['Lat_to'])
                single_event_as_region['longitude'] = (single_event['Lon_from'],\
                                                      single_event['Lon_to'])
                single_event_as_region['time'] = (single_event['Time_start'].to_pydatetime(),\
                                                      single_event['Time_stop'].to_pydatetime())
                self.__extremes_regions__.update({extreme_event_id : single_event_as_region})

            except KeyError:
                self.__logger__.error("Entry not found in catalogue. Please check spelling of input. These are the available entries: \n{0}".format('\n'.join(list(ex_table.index.values))))
                raise
                
        # we assume that the read in regions fully cover the event in space and time (TODO: build in a check for this)
        # TODO: still needed?
                
        self.__logger__.info(self.__extremes_regions__)
        
#        self.__extremes_regions__ = self.__regions__
        # Loop over the different regions (i.e. the different events)
        for r,def_r in self.__extremes_regions__.items():
            # Now define the three cubes. Note that now they are really cubes,
            # not dictionaries that contain cubes.

            # The event cube spans the event in space and time
            event_cube = self.__spatiotemp_subsets__(cube,{r:def_r})[r]

            # The clim cube spans the event in space, with the time
            # spanning all available timesteps in the dataset
            spatial_def = def_r.copy()
            spatial_def.pop("time")
            clim_cube = self.__spatiotemp_subsets__(cube,{r:spatial_def})[r]

            # The ex cube is created here, and used later for saving the extreme
            # climatology, it has the same shape as the event_cube
            ex_cube = event_cube.copy()
            
            # Now loop over each gridpoint in the selected region
            counter_gridpoints = 1
            n_gridpoints = event_cube.shape[1]*event_cube.shape[2]

            if num_processors>1:
                self.__logger__.info("Start multi process calculation of extreme climatology " +
                                     "for %s gridpoints using multiprocessing",n_gridpoints)
                # setting up a pool
                # TODO make sure that this is machine compatiple (or leave it to the user)
                pool = Pool(processes = num_processors)
                
                # get an iterator for the the positions
                positions = list(it.product(range(event_cube.shape[1]), range(event_cube.shape[2])))
                # Now loop in a processor distributed manner over the positions in the data
                # returns ex_list as percentile values with the number in the timeseries
                ex_list = pool.starmap(extremes_1D, 
                                            zip(positions,
                                                it.repeat(event_cube),
                                                it.repeat(clim_cube),
                                                it.repeat(ex_cube),
                                                it.repeat(window_size),
                                                it.repeat(which_percentile),
                                                it.repeat(min_measurements),
                                                )
                                        )
                # map the positions and the results back to the cube
                extremes_1d_redistribute(ex_cube,
                                         ex_list,
                                         positions)
                pool.close()
            else:
                self.__logger__.info("Start single process calculation of extreme climatology " +
                                     "for %s gridpoints without multiprocessing",n_gridpoints)
                # Now loop over the data
                for ii in range(event_cube.shape[1]):
                    for jj in range(event_cube.shape[2]):
                        # Convert this gridpoint to a pandas timeseries object
                        gridpoint_ts = ipd.as_series(clim_cube[:,ii,jj])
                        for n,doy in enumerate(event_cube.coord('day_of_year').points):
                            tmin = (doy - window_size) % 366
                            tmax = (doy + window_size) % 366
                            if not tmin:
                                tmin = 366
                            if not tmax:
                                tmax = 366
                            if tmin <= tmax:
                                doy_window = (tmin <= gridpoint_ts.index.dayofyear) &\
                                             (gridpoint_ts.index.dayofyear <= tmax)
                            else:
                                doy_window = (tmin <= gridpoint_ts.index.dayofyear) |\
                                             (gridpoint_ts.index.dayofyear <= tmax)
                            # Extract the right data points
                            gridpoint_sample = gridpoint_ts[doy_window]
                            # Check if there are enough valid measurements in the sample
                            if np.isfinite(gridpoint_sample).sum() > min_measurements:
                                perc_val = np.nanpercentile(gridpoint_ts[doy_window],which_percentile)
                            else:
                                perc_val = np.nan
                            ex_cube.data[n,ii,jj] = perc_val
                        self.__logger__.info("Progress of xclim: %s percent",\
                                             np.round(100.*(counter_gridpoints/n_gridpoints),decimals=1))
                        counter_gridpoints += 1

            #TODO think about propagation of nan values
            ex_cube.data = np.ma.masked_equal(ex_cube.core_data(), np.nan)
            
            self.__logger__.info("Finished calculation of extreme climatology")

            # Calculate exceedance of the extremes climatology
            if which_percentile > 50: # we are interested in over threshold values
                amplitude = event_cube - ex_cube
            elif which_percentile < 50: # we are interested in under threshold values
                amplitude = ex_cube - event_cube
            else:
                self.__loger__.error("Percentile can not be 50, that wouldn't be an extreme climatology")
                raise ValueError
            # Note that due to the above check, amplitude values of the event are always positive
            # therefore mask negative values
            amplitude.data = np.ma.masked_where(amplitude.core_data() <= 0, amplitude.core_data(), copy = True)
            # Set long name
            amplitude.long_name = 'Amplitude of '+event_cube.long_name
            
            # This check assures that the dimensions differ in size, otherwise
            # broadcasting in np.atleast_3d could fail and produce erroneous 
            # results without notice
            assert(sorted(list(set(amplitude.shape)))==sorted(list(amplitude.shape)))


            #####################################
            ###  calculate the three metrics  ###
            #####################################
            # severity
            severity = amplitude * np.atleast_3d(np.array([np.diff(bds) for bds in event_cube.coord("time").bounds]))
            severity.data = np.ma.masked_where(~np.isfinite(severity.core_data()),\
                                               severity.core_data(), copy = True)
            severity.units = cf_units.Unit(str(cube.units) + 
                                           " " + 
                                           str(cube.coord("time").units).split(" ")[0])
            severity = severity.collapsed("day_of_year", iris.analysis.SUM) 
            severity.long_name = "severity"

            # magnitude
            magnitude = amplitude.collapsed("day_of_year", iris.analysis.MAX)
            magnitude.data = np.ma.masked_where(~np.isfinite(magnitude.core_data()),\
                                               magnitude.core_data(), copy = True)
            magnitude.long_name = "magnitude"

            # duration
            duration = (amplitude * 0 + 1) * np.atleast_3d(np.array([np.diff(bds) for bds in event_cube.coord("time").bounds]))
            duration.data = np.ma.masked_where(~np.isfinite(duration.core_data()),\
                                               duration.core_data(), copy = True)
            duration.units = (str(cube.coord("time").units).split(" ")[0])
            duration = duration.collapsed("day_of_year", iris.analysis.SUM)
            duration.long_name = "duration"

            # Now calculate the spatial event mask from the severity as outlined in 
            # section 3.2 of the Extreme Catalogue C3S_D511.1.5_Extreme_Catalogue_V1.pdf
            event_mask_threshold = float(severity.collapsed(["latitude", "longitude"],
                                         iris.analysis.MEDIAN).data)
            event_mask = severity.data > event_mask_threshold
            # Here we move from a masked array to one single bool array 
            # where True values occur if threshold exceeded and input data (severity) 
            # was not masked.
            event_mask2d = event_mask.data & ~event_mask.mask

            # Now expand it over time
            event_mask3d = np.broadcast_to(event_mask2d, event_cube.shape)

            # Now create copies of the different event metrics and mask them with the event mask
            severity_masked = severity.copy()
            magnitude_masked = magnitude.copy()
            duration_masked = duration.copy()

            # Set the event_mask as mask for each of the three above
            severity_masked.mask = event_mask2d
            magnitude_masked.mask = event_mask2d
            duration_masked.mask = event_mask2d

            # calculate spatial averages 
            grid_areas = iris.analysis.cartography.area_weights(severity_masked)
            severity_av = severity_masked.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas)             
            magnitude_av = magnitude_masked.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas) 
            duration_av = duration_masked.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas)

            # now calculate time evolution of amplitude within defined event_mask
            amplitude_time = amplitude.copy()
            amplitude_time.mask = event_mask3d

            # Note that collapsing needs to be done in two steps, otherwise masked values
            # are propagated over the full array, which is not the preferred behaviour.
            amplitude_time = amplitude_time.collapsed(["latitude"],\
                                         iris.analysis.MEAN)
            amplitude_time = amplitude_time.collapsed(["longitude"],\
                                         iris.analysis.MEAN)

            # Plot time evolution
            plt.clf()
            qplt.plot(amplitude_time)
            plt.xticks(rotation=45)
            plt.savefig(self.__plot_dir__ + os.sep + r + "_amplitude_time_" + ".png")


            # calculate extent
            extent = ((duration * 0 + 1.) * grid_areas).collapsed(["latitude","longitude"], iris.analysis.SUM)/1e6
            extent.units = cf_units.Unit("km2")
            
            # now calculate time evolution of amplitude within defined event_mask
            extent_time = amplitude.copy()
            extent_time = (extent_time * 0 + 1.) * np.broadcast_to(grid_areas, extent_time.shape)/1e6
            extent_time.units = cf_units.Unit("km2")
            extent_time.long_name = "Extent (based on grid resolution)"
            extent_time.mask = event_mask3d

            # Note that collapsing needs to be done in two steps, otherwise masked values
            # are propagated over the full array, which is not the preferred behaviour.
            extent_time = extent_time.collapsed(["latitude"],\
                                         iris.analysis.MEAN)
            extent_time = extent_time.collapsed(["longitude"],\
                                         iris.analysis.SUM)

            # Plot time evolution
            plt.clf()
            qplt.plot(extent_time)
            plt.xticks(rotation=45)
            plt.savefig(self.__plot_dir__ + os.sep + r + "_extent_time_" + ".png")
            
            # set up table
            # TODO export csv
            self.__logger__.info("Extremes table")
            self.__logger__.info("mean severity: {:.2f} {}".format(severity_av.data, severity_av.units))
            self.__logger__.info("mean magnitude: {:.2f} {}".format(magnitude_av.data, magnitude_av.units))
            self.__logger__.info("mean duration: {:.2f} {}".format(duration_av.data, duration_av.units))
            self.__logger__.info("extent: {:.2f} {}".format(extent.data, extent.units))
            # plotting for trials
            
            for dat in ["severity", "magnitude", "duration"]:
                #TODO add event_mask_2d to the plots. 

                #self.__logger__.info(locals()[dat])
            
                qplt.pcolormesh(locals()[dat])#, vmin = 250, vmax = 300)
                
                # Add coastlines to the map created by contourf.
                plt.gca().coastlines("10m")
                plt.colorbar()
                plt.savefig(self.__plot_dir__ + os.sep + r + '_' + dat + ".png")
                
                plt.close()


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
    



def extremes_1D(ind, event_cube,  clim_cube, ex_cube, window_size, which_percentile, min_measurements):
    ii = ind[0]
    jj = ind[1]
    # Convert this gridpoint to a pandas timeseries object
    gridpoint_ts = ipd.as_series(clim_cube[:,ii,jj])
    # Get the doys that need to be processed
    perc_val_ts = []
    for n,doy in enumerate(event_cube.coord('day_of_year').points):
        tmin = (doy - window_size) % 366
        tmax = (doy + window_size) % 366
        if not tmin:
            tmin = 366
        if not tmax:
            tmax = 366
        if tmin <= tmax:
            doy_window = (tmin <= gridpoint_ts.index.dayofyear) &\
                         (gridpoint_ts.index.dayofyear <= tmax)
        else:
            doy_window = (tmin <= gridpoint_ts.index.dayofyear) |\
                         (gridpoint_ts.index.dayofyear <= tmax)
        # Extract the right data points
        gridpoint_sample = gridpoint_ts[doy_window]
        # Check if there are enough valid measurements in the sample
        if np.isfinite(gridpoint_sample).sum() > min_measurements:
            perc_val = np.nanpercentile(gridpoint_ts[doy_window],which_percentile)
        else:
            perc_val = np.nan
        
        # convert the individual values to a list of tuples 
        # providing percentiles and number in the time series
        perc_val_ts.append((perc_val,n))
        
    return perc_val_ts

def extremes_1d_redistribute(ex_cube, extrems_1D, positions):
    # puts back the data into the specific coordinates
    for ind,(ii,jj) in enumerate(positions):
        for perc_val, n in extrems_1D[ind]:
            ex_cube.data[n,ii,jj] = perc_val
    return

    
    
