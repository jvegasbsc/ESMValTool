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
import xarray
import cartopy.crs as ccrs

from .c3s_511_basic import Basic_Diagnostic_SP
from .libs.MD_old.ESMValMD import ESMValMD
from .libs.trends.trends_core3d import linear_trend, theilsen_trend, mannkendall
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
        
        list_of_plots = []
        
#        for rcoord in ["day_of_month", "day_of_year", "month_number", "year"]:
#            if rcoord in [coord.name() for coord in cube.coords()]:
#                try:
#                    cube.remove_coord(rcoord)
#                except:
#                    cube.remove_aux_factory(rcoord)
        
        # basic calculations
        lintrend,linpvalue = linear_trend(xarray.DataArray.from_iris(cube))
        theilsen_slope = theilsen_trend(xarray.DataArray.from_iris(cube))
        mk_result = mannkendall(xarray.DataArray.from_iris(cube))
        
        self.__logger__.info(lintrend)
        self.__logger__.info(linpvalue)
        self.__logger__.info(theilsen_slope)
        self.__logger__.info(mk_result)
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        lintrend.plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "lintrend.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        linpvalue.plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "pval.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        theilsen_slope.plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "theilsen.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        ((theilsen_slope-lintrend)/theilsen_slope).plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "trenddiff.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        mk_result.plot(ax=ax, transform=ccrs.PlateCarree(),robust=True, vmin=-1, vmax=1) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "mankendall.png")
        
        
        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_trends_input",
                                      expfile="_trends.txt",
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