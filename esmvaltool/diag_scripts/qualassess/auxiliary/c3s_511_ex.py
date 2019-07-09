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

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.predef.ecv_lookup_table import ecv_lookup
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
        if self.var3D is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(max(self.levels))
            dataset_id = [self.__dataset_id__[0], "at", "level", str(
                max(self.levels)), str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        ##################
        
        # we assume the regions to be spatio-temporally buffered.
        list_of_plots = []
        
        for r in self.__regions__:
            spat_r = self.__regions__[r]
            loc_cube = self.__spatiotemp_subsets__(cube,{r:spat_r})
            
            list_of_cubes = []
            
            for (entry,cube) in loc_cube.items():
                for rcoord in ["day_of_month", "month_number", "year"]:
                    if rcoord in [coord.name() for coord in cube.coords()]:
                        cube.remove_coord(rcoord)
            
            for yx_slice in loc_cube[r].slices(['time']):
                agg_cube = yx_slice.aggregated_by("day_of_year",iris.analysis.MEAN) 
                list_of_sub_cubes=[]
                for doy in np.sort(agg_cube.coords("day_of_year")[0].points):
                    loc_slice = agg_cube.extract(iris.Constraint(day_of_year=doy))
                    tmin = (doy - window_size) % 366
                    tmax = (doy + window_size) % 366
                    if not tmin:
                        tmin = 366
                    if not tmax:
                        tmax = 366
                    if tmin > tmax:
                        doy_sel = yx_slice.extract(
                                iris.Constraint(coord_values={'day_of_year':
                                    lambda cell: 1 <= cell <= tmin or
                                                 tmax <= cell <= 366}))
                    else:
                        doy_sel = yx_slice.extract(
                                iris.Constraint(coord_values={'day_of_year':
                                    lambda cell: tmin <= cell <= tmax}))
                    try:
                        perc = doy_sel.collapsed("time",iris.analysis.PERCENTILE,percent=[which_percentile]).core_data()
                    except:
                        perc = doy_sel.core_data()
                    if np.ma.is_masked(perc):
                        masked_val = perc.data
                        perc = perc.data
                            
                    loc_slice.data = perc
                    loc_slice.remove_coord("day_of_year")
                    
                    list_of_sub_cubes.append(loc_slice)
                    
                t_cube = iris.cube.CubeList(list_of_sub_cubes).merge()[0]
                    
                iris.util.promote_aux_coord_to_dim_coord(t_cube, "time")
                iris.util.new_axis(t_cube, "latitude")
                iris.util.new_axis(t_cube, "longitude")
                list_of_cubes.append(t_cube)
                
                
            clim_cube = iris.cube.CubeList(list_of_cubes).merge()[0]
            
            clim_cube.data = np.ma.masked_equal(clim_cube.core_data(), masked_val)
            
            self.__logger__.info(clim_cube)
            
            for yx_slice in clim_cube.slices(['latitude', 'longitude']):
                import iris.quickplot as qplt
                
                # Draw the contour with 25 levels.
                qplt.pcolormesh(yx_slice, vmin=240, vmax=300)
                
                # Add coastlines to the map created by contourf.
                plt.gca().coastlines()
                
                plt.savefig(self.__plot_dir__ + os.sep + str(yx_slice.coord("time").points[0]).replace(".","_") + ".png")
                
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