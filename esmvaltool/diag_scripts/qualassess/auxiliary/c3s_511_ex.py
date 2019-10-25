#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:59:59 2018

@authors: bmueller, bcrezee

Caveat: 

  - does not work yet on events that are symmetric in lat/lon (can not safely broadcast ERROR will be raised)

"""

import pandas as pd
import iris
import iris.pandas as ipd
import iris.plot as iplt
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime
import numpy as np
import cf_units
from sharedutils import parallel_apply_along_axis

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
from .libs.c3s_511_util import read_extreme_event_catalogue
from .plots.basicplot import Plot2D, Plot2D_blank, plot_setup
    
from multiprocessing import Pool
import itertools as it
import esmvalcore.preprocessor as pp

def get_doy_mask(cube, doy_range, window_size):
    if window_size%2==0:
        raise ValueError("Window size should be an uneven number")
    sample = cube[:,0,0]
    try:
        iris.coord_categorisation.add_day_of_year(sample, 'time', name='day_of_year')
    except ValueError:
        pass
    doys = sample.coord('day_of_year').points
    
    half_window = (window_size - 1) / 2
    
    doymask = np.empty((len(doy_range),)+doys.shape,dtype=bool)
    
    for n,doy in enumerate(doy_range):
        tmin = (doy - half_window) % 366
        tmax = (doy + half_window) % 366
        if not tmin:
            tmin = 366
        if not tmax:
            tmax = 366
        if tmin <= tmax:
            doy_window = (tmin <= doys) &\
                         (doys <= tmax)
        else:
            doy_window = (tmin <= doys) |\
                         (doys <= tmax)
        doymask[n,:] = doy_window
    return doymask


def get_xclim(cube, percentile, doy_range, window_size=None):
    print("Getting mask for day of year")
    doymask = get_doy_mask(cube, doy_range, window_size) # shape (len(doy_range), ntimesteps)
    result = np.empty(doymask.shape[:1] + cube.shape[1:])
    print(result.shape)
    for n,doy in enumerate(doy_range):
        print("Reading data into memory for window around day-of-year {0}".format(doy))
        doy_data = cube[doymask[n, :], :, :].data.filled(np.nan)
        print("Calculating clim...")
        result[n,:,:] = parallel_apply_along_axis(np.nanpercentile, 0, doy_data, percentile)
    return result


def subtract_xclim(cube, percentile=None, refcube=None):
    if not refcube:
        cube = refcube
    try:
        iris.coord_categorisation.add_day_of_year(cube, 'time', name='day_of_year')
    except ValueError:
        pass
    xclim_array = get_xclim(cube, percentile, cube.coord('day_of_year').points)
    # Now repeat
    if percentile >= 50:
        cube.data = cube.data - np.ma.masked_invalid(xclim_array)
    elif percentile < 50:
        cube.data = np.ma.masked_invalid(xclim_array) - cube.data
    return cube


class ex_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(ex_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
        
        # all required input can be extracted from the extremes dictionary
#        self.__logger__.info(self.__extremes__)
        
#        self.__regions__.update({
#            'Europe_1999': {
#                'latitude': (30, 75),
#                'longitude': (-10, 35),
#                'time': (datetime.datetime(1999, 5, 1),
#                         datetime.datetime(1999, 9, 30)
#                         )}})  # default region
        
        # Initialize extremes regions as empty
        self.__extremes_regions__ = dict()

    def _extremes_preprocessing(self, cube):
        regridres = '1x1'
        self.__logger__.info('Calculating daily means')
        #cube = pp.daily_statistics(cube)
        self.__logger__.info('Regridding to %s',regridres)
        #cube = pp.regrid(cube,regridres, scheme='area_weighted')
        return cube


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
        # self.__extremes__ provides the information given in the recipe as a dictionary
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
                # self.__extremes_regions__ provides the regions for the events given in the recipe as a dictionary (read from respective catalogue)
                self.__extremes_regions__.update({extreme_event_id : single_event_as_region})
            except KeyError:
                self.__logger__.error("Entry not found in catalogue. Please check spelling of input. These are the available entries: \n{0}".format('\n'.join(list(ex_table.index.values))))
                raise

        # Initialize a pandas dataframe for saving the table of metrics
        df_metrics = pd.DataFrame(columns=['severity','magnitude','duration','extent'],index=self.__regions__.keys(),dtype=float)

        # initialize a dict for the multiple data sets from regions
        extremes_measures = {}

        # Loop over the different regions (i.e. the different events)
        for r,def_r in self.__extremes_regions__.items():
            self.__logger__.info("Handling event {}".format(r))
            # Now define the three cubes. Note that now they are really cubes,
            # not dictionaries that contain cubes.

            # The event cube spans the event in space and time
            try:
                event_cube = self.__spatiotemp_subsets__(cube,{r:def_r})[r]
            except ValueError:
                self.__logger__.warning("The following extreme event was not \
                                        found in this dataset: {0}".format(r))
                continue

            # The clim cube spans the event in space, with the time
            # spanning all available timesteps in the dataset
            spatial_def = def_r.copy()
            spatial_def.pop("time")
            clim_cube = self.__spatiotemp_subsets__(cube,{r:spatial_def})[r]

            # Now apply common preprocessing steps to get to daily data at 1x1 degree resolution
            clim_cube = self._extremes_preprocessing(clim_cube)
            event_cube = self._extremes_preprocessing(event_cube)

            # The ex cube is created here, and used later for saving the extreme
            # climatology, it has the same shape as the event_cube
            ex_cube = event_cube.copy()
            
            # Now loop over each gridpoint in the selected region
            counter_gridpoints = 1
            n_gridpoints = event_cube.shape[1]*event_cube.shape[2]
            clim_cube.data
            self.__logger__.info("Data has been read into memory")
            self.__logger__.info("Start multi process calculation of extreme climatology " +
                                     "for %s gridpoints using multiprocessing",n_gridpoints)
            # Now calculate the xclim from clim_cube and subtract it from the event_cube.
            amplitude = subtract_xclim(event_cube, percentile=which_percentile, refcube=clim_cube)
            self.__logger__.info("Finished calculating amplitude")

            # Note that due to the above check, amplitude values of the event are always positive
            # therefore mask negative values
            amplitude.data = np.ma.masked_where(amplitude.core_data() <= 0, amplitude.core_data(), copy = True)
            # Set long name
            amplitude.long_name = 'Amplitude of '+event_cube.long_name
            
            # The below check assures that the dimensions differ in size, otherwise
            # broadcasting in np.atleast_3d could fail and produce erroneous 
            # results without notice
            #TODO fix the below
            try:
                assert(sorted(list(set(amplitude.shape)))==sorted(list(amplitude.shape)))
            except AssertionError:
                self.__logger__.error("Can not safely broadcast, since not all dimensions \
                                       of cube differ in size. Therefore skipping this event.")
                continue

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
            
            # log all results in a dictionary
            extremes_measures.update({r:
                {"severity":severity,
                 "magnitude":magnitude,
                 "duration":duration,
                 "time_series":{"amplitude":amplitude_time,
                                "extent":extent_time}}
                    }
                    )
            
            # Add metrics to pd dataframe
            df_metrics.loc[r] = pd.Series({'severity': severity_av.data, 'magnitude':magnitude_av.data, 'duration': duration_av.data, 'extent': extent.data})

        self.__logger__.info(extremes_measures)
        
        # calculate the common minimum and maximum for the plots
        common_minmax = {}
        
        for dat in ["severity", "magnitude", "duration"]:
            loc_list = []
            for reg, m in extremes_measures.items():
                loc_list.append(np.nanpercentile(m[dat].data.compressed(),[5, 95]))
            allmin = np.nanmin([np.nanmin(loc_el) for loc_el in loc_list])
            allmax = np.nanmax([np.nanmax(loc_el) for loc_el in loc_list])
            common_minmax.update({dat:[allmin,allmax]})
            
        for r,m in extremes_measures.items():
            # plotting maps
            for dat in ["severity", "magnitude", "duration"]:
                #TODO add event_mask_2d to the plots.
                
                filename = self.__plot_dir__ + os.sep + \
                    basic_filename + \
                    "_" + r + "_" + dat + \
                    "." + self.__output_type__
                list_of_plots.append(filename)
    
                try:
                    x = Plot2D(m[dat])
    
                    caption = str(dat.title() +
                                  ' maps of ' +
                                  ecv_lookup(self.__varname__) +
                                  ' for the extreme event ' + r +
                                  ' for the data set ' +
                                  " ".join(dataset_id) + ' (' +
                                  self.__time_period__ + ').')
    
                    fig = plt.figure()
    
                    (fig, ax, caption) = plot_setup(d="time",
                                                    numfigs=1,
                                                    fig=fig,
                                                    caption=caption)

                    vminmax = common_minmax[dat]
    
                    x.plot(ax=ax,
                           color={"Extremes": "YlOrRd"},
                           color_type="Extremes",
                           title=" ".join([self.__dataset_id__[indx] for
                                           indx in [0, 2, 1, 3]]) + \
                                 " (" + self.__time_period__ + ")",
                           vminmax=vminmax,
                           ext_cmap="both",
                           dat_log = self.log_data)
                    fig.savefig(filename)
                    plt.close(fig.number)
    
                    ESMValMD("meta",
                             filename,
                             self.__basetags__ + 
                             ['DM_regional', 'C3S_extremes'],
                             caption +  
                             ('NA-values and values of 0 and below are shown ' + 
                              'in grey.' if self.log_data else
                              'NA-values are shown in grey.'),
                             '#C3S' + "extremes" + dat + \
                             self.__varname__,
                             self.__infile__,
                             self.diagname,
                             self.authors)
    
                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(
                        exc_tb.tb_frame.f_code.co_filename)[1]
                    self.__logger__.error(
                        exc_type, fname, exc_tb.tb_lineno)
                    self.__logger__.error(dat)
                    self.__logger__.error('Warning: blank figure!')
    
                    x = Plot2D_blank(m[dat])
    
                    fig = plt.figure()
    
                    (fig, ax, caption) = plot_setup(d="time", 
                                                    numfigs=1,
                                                    fig=fig,
                                                    caption=caption)
    
                    x.plot(ax=ax,
                           color={"Extremes": "YlOrRd"},
                           color_type="Extremes",
                           title=" ".join([self.__dataset_id__[indx] for 
                                           indx in [0, 2, 1, 3]]) + \
                                 " (" + self.__time_period__ + ")")
                    fig.savefig(filename)
                    plt.close(fig.number)
    
                    ESMValMD("meta",
                             filename,
                             self.__basetags__ +
                             ['DM_global', 'C3S_extremes'],
                             caption +
                             '; Data can ' +
                             'not be displayed due to cartopy error!',
                             '#C3S' + "extrenes" + dat +
                             self.__varname__,
                             self.__infile__, 
                             self.diagname, 
                             self.authors)
                    
            for dat in ["time_series"]: 
                amplitude_time = m[dat]["amplitude"]
                extent_time = m[dat]["extent"]
                if all(amplitude_time.data.mask) or all(extent_time.data.mask):
                    self.__logger__.warning("Extremes: spatially aggregated amplitude or extent are all masked! No lineplots produced.")
                else:
                    # plotting lineplots
                    filename = self.__plot_dir__ + os.sep + \
                        basic_filename + \
                        "_" + r + "_" + "temp_extent_amplitude" + \
                        "." + self.__output_type__
                    list_of_plots.append(filename)
        
                    caption = str('Temporal progress of amplitude and extent of ' +
                                  ecv_lookup(self.__varname__) +
                                  ' for the extreme event ' + r +
                                  ' for the data set ' +
                                  " ".join(dataset_id) + ' (' +
                                  self.__time_period__ + ').')
        
                    fig = plt.figure()
                    fig.set_figwidth(1.7 * fig.get_figwidth())
                    fig.set_figheight(2.2 * fig.get_figheight())
        
                    gs = gridspec.GridSpec(8, 1)
                    ax = np.array([plt.subplot(gs[0:4,:]),plt.subplot(gs[4:8,:])])
                    plt.sca(ax[0])
                    iplt.plot(amplitude_time)
                    plt.ylabel(amplitude_time.long_name + " [{}]".format(amplitude_time.units))
                    plt.tick_params(axis='x',
                                    which='both',      
                                    bottom=False,      
                                    top=True,         
                                    labelbottom=False) 
                    plt.grid(color="k", linestyle=':')
                    plt.sca(ax[1])
                    iplt.plot(extent_time)
                    plt.ylabel(extent_time.long_name + " [{}]".format(extent_time.units))
                    plt.xticks(rotation=45)
                    plt.grid(color="k", linestyle=':')
                    fig.align_ylabels()
                    plt.tight_layout()
        
                    fig.savefig(filename)
                    plt.close(fig.number)
        
                    ESMValMD("meta",
                             filename,
                             self.__basetags__ + 
                             ['DM_regional', 'C3S_extremes'],
                             caption,
                             '#C3S' + "extremes" + "ExtAmp" + \
                             self.__varname__,
                             self.__infile__,
                             self.diagname,
                             self.authors)

        # Add units to column names
        column_rename_dict = {'severity': 'severity [{0}]'.format(severity_av.units), 'magnitude': 'magnitude [{0}]'.format(magnitude_av.units), 'duration': 'duration [{0}]'.format(duration_av.units), 'extent': 'extent [{0}]'.format(extent.units)}
        df_metrics = df_metrics.rename(columns=column_rename_dict)
        df_metrics = df_metrics.astype(float)

        # Saving of the table to csv and html (saved to work. plot_dir only contains plots)}
        savename = self.__work_dir__ + os.sep + basic_filename + "_extreme_event_metrics."
        savename_csv = savename + "csv"
        savename_html = savename + "html"

        # Save to csv (keeping full precision)
        self.__logger__.info("Saving metric csv-table as: {0}".format(savename_csv))
        df_metrics.to_csv(savename_csv)

        # Save to html with specified precision
        html_metrics = df_metrics.style.set_precision(2).render()
        self.__logger__.info("Saving metric html-table as: {0}".format(savename_html))
        with open(savename_html,mode='w+') as handle:
            handle.write(html_metrics)
            
        # end multiprocessing
        if num_processors>1:
        # closing the pool
            pool.close()

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
