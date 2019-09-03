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
from .libs.trend_framework.trends_core3d import linear_trend, theilsen_trend, mannkendall
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
        
        xcube = xarray.DataArray.from_iris(cube)
        
        # basic calculations
        lintrend,linpvalue = linear_trend(xcube)
        theilsen_slope = theilsen_trend(xcube)
        mk_result = mannkendall(xcube)
        
        # prepare cube conversion linear trend
        if (" " in lintrend.name):
            lintrend.name = lintrend.name.strip(" ")[0]
        lintrend.attrs['units'] = lintrend.attrs['units'].replace("timestep", "({} {})".format(self.__avg_timestep__[1],cube.coord("time").units.name.split(" ")[0]))
        lintrend.attrs.update({'standard_name': None})
        lintrend.attrs.update({'long_name': "Linear Trend of {}".format(xcube.attrs["long_name"])})
        lintrend.attrs.update({'cell_methods': 'time: linear trend'})
        
        # prepare cube conversion pvalue
        if (" " in linpvalue.name):
            linpvalue.name = linpvalue.name.strip(" ")[0]
        linpvalue.attrs.update({'standard_name': None})
        linpvalue.attrs.update({'long_name': "Linear Trend of {} (p-value)".format(xcube.attrs["long_name"])})
        linpvalue.attrs.update({'cell_methods': 'time: linear trend (pvalue)'})
        
        # prepare cube conversion theilsen
        if (" " in theilsen_slope.name):
            theilsen_slope.name = theilsen_slope.name.strip(" ")[0]
        theilsen_slope.attrs['units'] = theilsen_slope.attrs['units'].replace("per timestep", "/ ({} {})".format(self.__avg_timestep__[1],cube.coord("time").units.name.split(" ")[0]))
        theilsen_slope.attrs.update({'standard_name': None})
        theilsen_slope.attrs.update({'long_name': "Theil-Sen Trend of {}".format(xcube.attrs["long_name"])})
        theilsen_slope.attrs.update({'cell_methods': 'time: Theil-Sen trend (pvalue)'})
        
        # prepare cube conversion mankendall
        if (" " in mk_result.name):
            mk_result.name = mk_result.name.strip(" ")[0]
        mk_result.attrs.update({'standard_name': None})
        mk_result.attrs.update({'long_name': "Sign of Significant Trend of {}".format(xcube.attrs["long_name"])})
        mk_result.attrs.update({'cell_methods': 'time: mankendall'})
        
        # adjust attributes
        for key,val in xcube.attrs.items():
            if key not in lintrend.attrs.keys():
                lintrend.attrs.update({key: val})
            if key not in linpvalue.attrs.keys():
                linpvalue.attrs.update({key: val})
            if key not in theilsen_slope.attrs.keys():
                theilsen_slope.attrs.update({key: val})
            if key not in mk_result.attrs.keys():
                mk_result.attrs.update({key: val})
                
        #conversions to cubes
        lintrend_cube = lintrend.to_iris()
        lintrend_cube.convert_units('{} / (10 years)'.format(cube.units))
        
        linpvalue_cube = linpvalue.to_iris()
        
        theilsen_slope_cube = theilsen_slope.to_iris()
        theilsen_slope_cube.convert_units('{} / (10 years)'.format(cube.units))
        
        mk_result_cube = mk_result.to_iris()
                
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        xarray.DataArray.from_iris(lintrend_cube).plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "lintrend.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        xarray.DataArray.from_iris(linpvalue_cube).plot(ax=ax, transform=ccrs.PlateCarree(),robust=True, vmin=0, vmax=1) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "pval.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        xarray.DataArray.from_iris(theilsen_slope_cube).plot(ax=ax, transform=ccrs.PlateCarree(),robust=True) # Note that computation takes place at this place
        ax.coastlines()
        plt.tight_layout()
        fig.savefig(self.__plot_dir__ + os.sep + "theilsen.png")
        
        fig = plt.figure(figsize=(15, 7))
        ax = fig.add_subplot(111, projection=ccrs.Robinson())
        xarray.DataArray.from_iris(mk_result_cube).plot(ax=ax, transform=ccrs.PlateCarree(),robust=True, vmin=-1, vmax=1) # Note that computation takes place at this place
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