#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 13:59:59 2018

@author: bmueller
"""

import iris
import os
import sys
import subprocess
import matplotlib.pyplot as plt
import datetime
import numpy as np
import xarray
import cartopy.crs as ccrs

LD_Lib_Path = os.environ.pop("LD_LIBRARY_PATH","")
with open(os.devnull) as devnull:
    which_R = subprocess.check_output(["which", "R"], stderr=devnull)
    which_R = which_R.decode("utf-8").replace("\n","/lib")
os.environ["LD_LIBRARY_PATH"] = "{}:{}".format(which_R,LD_Lib_Path)

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.trend_framework.c3s_511_trends import TrendLims1D
from .libs.predef.ecv_lookup_table import ecv_lookup
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup

class trnd_Diagnostic_SP(Basic_Diagnostic_SP):
    """
    class to implement additional diagnostics
    """
    
    def set_info(self, **kwargs):
        
        super(trnd_Diagnostic_SP, self).set_info(**kwargs)
        
        # add a region to the regions object
#        self.__regions__.update({
#            'MAR_region_Aug/Sept': {
#                'latitude': (-60, 50),
#                'longitude': (-60, 0),
#                'time': (datetime.datetime(2000, 8, 1),
#                         datetime.datetime(2000, 9, 30)
#                         )
#                }})
    
    def run_diagnostic(self):
#        self.sp_data = self.__spatiotemp_subsets__(self.sp_data)['Europe_2000']
        self.__do_trends__()

        self.__do_full_report__()

    def __do_trends__(self):
        
        this_function =  "trends example"
        
        # this the trends example
        
        self.__logger__.info(self.var3D)
        
        # handling of the cube restrictions copied from the basic diagnostic
        # needs to be set locally
        #################
        if not self.var3D:
            cube = self.sp_data
        else:
            cube = self.sp_data.extract(iris.Constraint(
                    coord_values={str(self.level_dim): lambda cell: cell == max(self.levels)}))

        # adjustment to ids and filenames
        if self.var3D not in [None, False]:
            basic_filename = self.__basic_filename__ + "_lev" + str(max(self.levels))
            dataset_id = [self.__dataset_id__[0], "at", "level", str(
                max(self.levels)), str(self.sp_data.coord(self.level_dim).units)]
    
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        ##################
        
        list_of_plots = []
        
#        for rcoord in ["day_of_month", "day_of_year", "month_number", "year"]:
#            if rcoord in [coord.name() for coord in cube.coords()]:
#                try:
#                    cube.remove_coord(rcoord)
#                except:
#                    cube.remove_aux_factory(rcoord)
        
#        xcube = xarray.DataArray.from_iris(cube)
        
        for ts_cube in cube.slices(["time"]):
            ts_cube = iris.util.new_axis(ts_cube, "latitude")
            ts_cube = iris.util.new_axis(ts_cube, "longitude")
            ts_cube.transpose([2,1,0])
            self.__logger__.info(ts_cube)
#            trend_obj = TrendLims1D("local")
#            trend_obj.initialize_through_realization_of_cube(ts_cube,0,0)
#            trend_obj.resample('Y')
#            trend_obj.do_trends()
#            self.__logger__.info(trend_obj.data_ts)
#            assert False, "development"
        
#        fig = plt.figure(figsize=(15, 7))
#        ax = fig.add_subplot(111, projection=ccrs.Robinson())
#        xarray.DataArray.from_iris(mk_result_cube).plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
#        ax.coastlines()
#        plt.tight_layout()
#        fig.savefig(self.__plot_dir__ + os.sep + "mankendall.png")
#        

        
        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_mytrends_input",
                                      expfile="_mytrends.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:

            self.reporting_structure.update(
                    {"trends": 
                        {"plots": list_of_plots,
                         "freetext": expected_input}})
        else:

            self.reporting_structure.update(
                    {"trends": 
                        {"plots": list_of_plots}})