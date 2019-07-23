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
        if self.var3D is not None:
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
            
            list_of_cubes = []
            
            # Slicing and merging cubes sometimes leads to the wrong dimensions
            # being propagated to coordinates instead of being auxiliary coordinates.
            # Therefore they are removed. 
            # TODO: handle this more specifically, in certain cases this might fail.
            for (entry,ecube) in loc_cube.items():
                for rcoord in ["day_of_month", "month_number", "year", self.level_dim]:
                    if rcoord in [coord.name() for coord in ecube.coords()]:
                        ecube.remove_coord(rcoord)
            
            # Save the bounds
            latbnds = ecube.coord("latitude").bounds
            lonbnds = ecube.coord("longitude").bounds
            
            # Now loop over each time step
            for yx_slice in loc_cube[r].slices(['time']):
                # Prepare for calculating extremes climatology
                agg_cube = yx_slice.aggregated_by("day_of_year",iris.analysis.MEAN) 
                list_of_sub_cubes=[]
                
                # Now loop over each doy and extract a temporal window surrounding this doy
                for doy in np.sort(agg_cube.coord("day_of_year").points):
                    # loc_slice is created here for carrying on the metadata
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
                    # BAS: this try except statement should be improved. 
                    # when does an exception occur, and why? these statements
                    # are terrible for debugging
                    if len(doy_sel.coord("time").points)>1:
                        perc = doy_sel.collapsed("time",
                                                 iris.analysis.PERCENTILE,
                                                 percent=[which_percentile]
                                                 ).core_data()
                    else:
                        perc = doy_sel.core_data()
                    # BAS: when does this happen?
                    if np.ma.is_masked(perc):
                        masked_val = perc.mask
                        perc = perc.data
                            
                    loc_slice.data = perc
                    loc_slice.remove_coord("day_of_year")
                    
                    list_of_sub_cubes.append(loc_slice)
                    
                # TODO: make make more explicit, merge('time')
                # t_cube can be interpreted as a timeseries of one pixel 
                t_cube = iris.util.squeeze(iris.cube.CubeList(list_of_sub_cubes).merge()[0])
                
                iris.util.promote_aux_coord_to_dim_coord(t_cube, "time")
                # Since there is no promote_scalar_coord_to_dim_coord, do it like this:
                iris.util.new_axis(t_cube, "latitude")
                iris.util.new_axis(t_cube, "longitude")
                list_of_cubes.append(t_cube)
                
                
            # BAS: Maybe assertions are not needed, since Iris handles data quite strict

            # Here Iris puts back the list of individual gridpoint timeseries 
            # on a 'lat lon grid'
            list_of_cubes_merged = iris.cube.CubeList(list_of_cubes).merge()
            assert(len(list_of_cubes_merged)==1)
            clim_cube = list_of_cubes_merged[0]
            
            clim_cube.data = np.ma.masked_equal(clim_cube.core_data(), masked_val)
            
            # This is needed to have the same order of dim_coords on both cubes
            cc_dim_order = [c.name() for c in clim_cube.coords()][:len(clim_cube.shape)]
            lc_dim_order = [c.name() for c in loc_cube[r].coords()][:len(loc_cube[r].shape)]
            clim_cube.transpose([cc_dim_order.index(i) for i in lc_dim_order])

            # BAS: check issue. At this point, loc_cube[r] still has aux coord day_of_year, even when removing it above
            # Here the actual exceedance of the extreme climatology is calculated
            incident_data = loc_cube[r]-clim_cube
            incident_data.data = np.ma.masked_where(incident_data.core_data() <= 0, incident_data.core_data(), copy = True)
            # Sidenote: 
            # TODO: Insert a check that incident_data cube does not have two or more dims with the same shape 
            # otherwise broadcast could fail and produce erroneous results without notice
            severity = incident_data * np.atleast_3d(np.array([np.diff(bds) for bds in incident_data.coord("time").bounds]))
            severity.units = cf_units.Unit(str(severity.units) + 
                                           " " + 
                                           str(severity.coord("time").units).split(" ")[0])
            severity = severity.collapsed("time", iris.analysis.SUM) 
            severity.long_name = "severity"
            
            magnitude = incident_data.collapsed("time", iris.analysis.MAX)
            magnitude.long_name = "magnitude"
            
            duration = (incident_data * 0 + 1) * np.atleast_3d(np.array([np.diff(bds) for bds in incident_data.coord("time").bounds]))
            duration.units = (str(duration.coord("time").units).split(" ")[0])
            duration = duration.collapsed("time", iris.analysis.SUM)
            duration.long_name = "duration"
            
            grid_areas = iris.analysis.cartography.area_weights(severity)
            severity_av = severity.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas)             
            magnitude_av = magnitude.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas) 
            duration_av = duration.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=grid_areas) 
            
            
            def ellipsoidal_distance(lat1, long1, lat2, long2):
                # https://www.johndcook.com/blog/2018/11/24/spheroid-distance/

                a = 6378137.0 # equatorial radius in meters 
                f = 1/298.257223563 # ellipsoid flattening 
                b = (1 - f)*a 
                tolerance = 1e-11 # to stop iteration
            
                phi1, phi2 = lat1, lat2
                U1 = arctan((1-f)*tan(phi1))
                U2 = arctan((1-f)*tan(phi2))
                L1, L2 = long1, long2
                L = L2 - L1
            
                lambda_old = L + 0
            
                while True:
                
                    t = (cos(U2)*sin(lambda_old))**2
                    t += (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda_old))**2
                    sin_sigma = t**0.5
                    cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda_old)
                    sigma = arctan2(sin_sigma, cos_sigma) 
                
                    sin_alpha = cos(U1)*cos(U2)*sin(lambda_old) / sin_sigma
                    cos_sq_alpha = 1 - sin_alpha**2
                    if cos_sq_alpha == 0: # for 0 there is no solution (=>nan)
                        cos_sq_alpha = tolerance/100
                    cos_2sigma_m = cos_sigma - 2*sin(U1)*sin(U2)/cos_sq_alpha
                    C = f*cos_sq_alpha*(4 + f*(4-3*cos_sq_alpha))/16
                
                    t = sigma + C*sin_sigma*(cos_2sigma_m + C*cos_sigma*(-1 + 2*cos_2sigma_m**2))
                    lambda_new = L + (1 - C)*f*sin_alpha*t
                    
                    if abs(lambda_new - lambda_old) <= tolerance:
                        break
                    else:
                        lambda_old = lambda_new
            
                u2 = cos_sq_alpha*((a**2 - b**2)/b**2)
                A = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320 - 175*u2)))
                B = (u2/1024)*(256 + u2*(-128 + u2*(74 - 47*u2)))
                t = cos_2sigma_m + 0.25*B*(cos_sigma*(-1 + 2*cos_2sigma_m**2))
                t -= (B/6)*cos_2sigma_m*(-3 + 4*sin_sigma**2)*(-3 + 4*cos_2sigma_m**2)
                delta_sigma = B * sin_sigma * t
                s = b*A*(sigma - delta_sigma)
            
                return s
            
            extent = (duration * 0 + 1.).collapsed("longitude", iris.analysis.SUM)
            
            pix_sizes_km2 = []
            for labs in latbnds:
                lobs = lonbnds[0] # nominally equal for all longitude bounds
                
                #approximated pixel side lengths (lat is the same for both sides)
                e1lon = ellipsoidal_distance(deg2rad(labs[1]), deg2rad(lobs[0]), deg2rad(labs[1]), deg2rad(lobs[1]))
                e2lon = ellipsoidal_distance(deg2rad(labs[0]), deg2rad(lobs[0]), deg2rad(labs[0]), deg2rad(lobs[1]))
                e0lat = ellipsoidal_distance(deg2rad(labs[0]), deg2rad(lobs[0]), deg2rad(labs[1]), deg2rad(lobs[0]))
                
                pix_size = round(e0lat * (e1lon + e2lon) / 2) #size rounded to m**2 assuming trapezoid
                
                pix_sizes_km2.append(pix_size/10e6)
            
            extent = (extent * np.array(pix_sizes_km2)).collapsed("latitude", iris.analysis.SUM)
            extent.units = cf_units.Unit("km2")
            
            extent2 = ((duration * 0 + 1.) * grid_areas).collapsed(["latitude","longitude"], iris.analysis.SUM)/10e6
            extent2.units = cf_units.Unit("km2")
            
            self.__logger__.info("Extremes table")
            self.__logger__.info("mean severity: {:.2f} {}".format(severity_av.data, severity_av.units))
            self.__logger__.info("mean magnitude: {:.2f} {}".format(magnitude_av.data, magnitude_av.units))
            self.__logger__.info("mean duration: {:.2f} {}".format(duration_av.data, duration_av.units))
            self.__logger__.info("extent: {:.2f} {}".format(extent.data, extent.units))
            self.__logger__.info("extent2: {:.2f} {}".format(extent2.data, extent2.units))
            
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
