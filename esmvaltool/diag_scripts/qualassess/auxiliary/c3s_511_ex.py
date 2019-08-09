#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:59:59 2018

@author: bmueller
"""

import iris
import os
import sys
import matplotlib.pyplot as plt
import datetime
import numpy as np
import cf_units
from scipy import sin, cos, tan, arctan, arctan2, arccos, pi, deg2rad

from multiprocessing import Pool
from itertools import repeat

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
#        self.__regions__ = dict({
#            'MAR_region_Jul/Sept': {
#                'latitude': (-60, 50),
#                'longitude': (-60, 0),
#                'time': (datetime.datetime(2000, 7, 1),
#                         datetime.datetime(2000, 9, 30)
#                         )
#                }})
    
    def run_diagnostic(self):
#        self.sp_data = self.__spatiotemp_subsets__(self.sp_data)['Europe_2000']
        self.__do_extremes__()

        self.__do_full_report__()

    def __do_extremes__(self):
        
        this_function =  "extremes example"
        
        which_percentile = 90
        window_size = 40 # one directional 5 => 11
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
            # Select this specific region, note that loc_cube is a dictionary with loc_cube[key]
            # holding the actual cubes.
            loc_cube = self.__spatiotemp_subsets__(cube,{r:def_r})

            spatial_def = def_r.copy()
            spatial_def.pop("time")
            spatial_loc_cube = self.__spatiotemp_subsets__(cube,{r:spatial_def})
            
            list_of_cubes = []
            
            #import IPython;IPython.embed()
            # Slicing and merging cubes sometimes leads to the wrong dimensions
            # being propagated to coordinates instead of being auxiliary coordinates.
            # Therefore they are removed. 
            # TODO: handle this more specifically, in certain cases this might fail.
            for (entry,ecube) in spatial_loc_cube.items():
                # Realize the data at this stage
                ecube.data
                for rcoord in ["day_of_month", "month_number", "year", self.level_dim]:
                    if rcoord in [coord.name() for coord in ecube.coords()]:
                        ecube.remove_coord(rcoord)
                        loc_cube[r].remove_coord(rcoord)
            

            n_gridpoints = spatial_loc_cube[r].coord('latitude').shape[0] *\
                spatial_loc_cube[r].coord('longitude').shape[0]
            self.__logger__.info("Start calculation of extreme climatology " +\
                                 "for %s gridpoints",n_gridpoints)
            
            # Now loop over each gridpoint in the selected region
            for nn,yx_slice in enumerate(spatial_loc_cube[r].slices(['time'])):
                # Prepare for calculating extremes climatology
                self.__logger__.info("Progress of xclim: %s percent",np.round(100.*((nn+1)/n_gridpoints),decimals=1))
                agg_cube = yx_slice.aggregated_by("day_of_year",iris.analysis.MEAN)
                list_of_sub_cubes=[]
                
                # Now loop over each doy and extract a temporal window surrounding this doy
                for doy in np.sort(agg_cube.coord("day_of_year").points):
                    # loc_slice is created here for carrying on the metadata
                    self.__logger__.info("#")
                    loc_slice = agg_cube.extract(iris.Constraint(day_of_year=doy))
                    tmin = (doy - window_size) % 366
                    tmax = (doy + window_size) % 366
                    self.__logger__.info([tmin,tmax])
                    if not tmin:
                        tmin = 366
                    if not tmax:
                        tmax = 366
                    if tmin > tmax:
                        doy_sel = yx_slice.extract(
                                iris.Constraint(coord_values={'day_of_year':
                                    lambda cell: 1 <= cell <= tmax or
                                                 tmin <= cell <= 366}))
                    else:
                        doy_sel = yx_slice.extract(
                                iris.Constraint(coord_values={'day_of_year':
                                    lambda cell: tmin <= cell <= tmax}))
    
                    self.__logger__.info(doy_sel)

                    # Do a check if there is actually data
                    if len(doy_sel.coord("time").points)>1:
                        self.__logger__.info(doy_sel.data)
                        self.__logger__.info(doy_sel.coord("time"))
                        perc = doy_sel.collapsed("time",
                                                 iris.analysis.PERCENTILE,
                                                 percent=[which_percentile]
                                                 ).core_data()
                    else:
                        perc = doy_sel.core_data()

                    # Catch the case of a masked array
                    if np.ma.is_masked(perc):
                        masked_val = perc.mask
                        perc = perc.data
                            
                    # Now the data is inserted back into a cube
                    loc_slice.data = perc
                    
                    list_of_sub_cubes.append(loc_slice)
                    
#                 TODO: make make more explicit, merge('day_of_year')
#                 t_cube can be interpreted as a timeseries of one pixel 

                t_cube = iris.util.squeeze(iris.cube.CubeList(list_of_sub_cubes).merge()[0])
                t_cube = cube_sorted(t_cube,"day_of_year")
                iris.util.promote_aux_coord_to_dim_coord(t_cube, "day_of_year")

                # Since there is no promote_scalar_coord_to_dim_coord, do it like this:
                iris.util.new_axis(t_cube, "latitude")
                iris.util.new_axis(t_cube, "longitude")
                list_of_cubes.append(t_cube)

            self.__logger__.info("Finished calculation of extreme climatology")

                
            # Here Iris puts back the list of individual gridpoint timeseries 
            # on a 'lat lon grid'
            self.__logger__.info("Start merging of gridpoints back to grid")
            list_of_cubes_merged = iris.cube.CubeList(list_of_cubes).merge()
            self.__logger__.info("Finished merging of gridpoints back to grid")
            assert(len(list_of_cubes_merged)==1)
            clim_cube = list_of_cubes_merged[0]
            
            clim_cube.data = np.ma.masked_equal(clim_cube.core_data(), masked_val)
            
            # This is needed to have the same order and names of dim_coords on both cubes
            self.__logger__.info("Reordering dimensions")
            cc_dim_order = [c.name() for c in clim_cube.coords()][:len(clim_cube.shape)]
            loc_cube_doy = loc_cube[r].copy()
            iris.util.promote_aux_coord_to_dim_coord(loc_cube_doy, "day_of_year")
            lc_dim_order = [c.name() for c in loc_cube_doy.coords()][:len(loc_cube_doy.shape)]
            clim_cube.transpose([cc_dim_order.index(i) for i in lc_dim_order])
            self.__logger__.info("Finished dimensions")

            # BAS: check issue. At this point, loc_cube[r] still has aux coord day_of_year, even when removing it above
            # Here the actual exceedance of the extreme climatology is calculated
            loc_cube_doy.remove_coord("time")
            clim_cube.remove_coord("time")
            doy_subset = [clim_cube.coord("day_of_year").points.tolist().index(i) for i in loc_cube_doy.coord("day_of_year").points.tolist()]
            
            incident_data = loc_cube_doy-clim_cube[doy_subset,:,:]
            incident_data.data = np.ma.masked_where(incident_data.core_data() <= 0, incident_data.core_data(), copy = True)
            
            # Sidenote: 
            # TODO: Insert a check that incident_data cube does not have two or more dims with the same shape 
            # otherwise broadcast could fail and produce erroneous results without notice
            
            # calculate severity
            severity = incident_data * np.atleast_3d(np.array([np.diff(bds) for bds in loc_cube[r].coord("time").bounds]))
            severity.units = cf_units.Unit(str(cube.units) + 
                                           " " + 
                                           str(cube.coord("time").units).split(" ")[0])
            severity = severity.collapsed("day_of_year", iris.analysis.SUM) 
            severity.long_name = "severity"
            
            # calculate magnitude
            magnitude = incident_data.collapsed("day_of_year", iris.analysis.MAX)
            magnitude.long_name = "magnitude"
            
            # calculate duration
            duration = (incident_data * 0 + 1) * np.atleast_3d(np.array([np.diff(bds) for bds in loc_cube[r].coord("time").bounds]))
            duration.units = (str(cube.coord("time").units).split(" ")[0])
            duration = duration.collapsed("day_of_year", iris.analysis.SUM)
            duration.long_name = "duration"
            
            # calculate spatial averages
            grid_areas = iris.analysis.cartography.area_weights(severity)
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
                self.__logger__.info(locals()[dat])
            
                qplt.pcolormesh(locals()[dat])#, vmin = 250, vmax = 300)
                
                # Add coastlines to the map created by contourf.
                plt.gca().coastlines()
                plt.colorbar()
                plt.savefig(self.__plot_dir__ + os.sep + dat + ".png")
                
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
