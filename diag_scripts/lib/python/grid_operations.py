# -*- coding: utf-8 -*-

"""
###############################################################################
GENERAL ROUTINES FOR OPERATIONS ON A (RECTILINEAR) GRID
###############################################################################
"""


from constants import Constant

import datetime
import netCDF4 as nc
import numpy as np


class GridOperations(object):
    """
    Description
        Class containing useful operations on a grid. Its instances are
        initialized by a single netCDF file which should contain temporal and
        spatial information of one variable.
        Please consider extending existing routines before adding new ones.
        Check the header of each method for documentation.

    Contents
        helper method reset
        helper method check_latlon_available
        helper method check_time_avaiable
        helper method select_region
        helper method gridcell_area
        helper method map_area
        helper method spatial_average
        helper method temporal_average
        method average
    """

    ###########################################################################

    def __init__(self, ncfile_path, variable):
        """
        Arguments
            ncfile_path : Path to the netCDF file
            variable    : Main variable which should be analyzed

        Description
            Initializes class instances.

        Modification history
            20171114-A_schl_ma: written
        """

        # Check if file exists
        try:
            self.ncfile = nc.Dataset(ncfile_path, "r")
        except:
            raise IOError("Invalid input: could not open the netCDF file " + \
                          " at {0}".format(ncfile_path))

        # Check if variable is valid
        try:
            valid_var = (type(variable) == str)
            self.variable = variable
        except:
            raise TypeError("Invalid input: variable")
        if (not valid_var):
            raise TypeError("Invalid input: variable has to be a string")

        # Check contents of netCDF file
        try:
            self.var_nc = self.ncfile[variable]
        except:
            raise ValueError("Invalid input: the given netCDF file does " + \
                             "not have a variable '{0}'".format(variable))
        try:
            self.time_nc = self.ncfile["time"]
            self.time_available = True
        except:
            self.time_available = False
        try:
            self.lat_nc = self.ncfile["lat"]
            self.lat_available = True
        except:
            self.lat_available = False
        try:
            self.lon_nc = self.ncfile["lon"]
            self.lon_available = True
        except:
            self.lon_available = False

    ###########################################################################

    def reset(self):
        """
        Arguments
            None

        Return value
            None

        Description
            Resets all internal arrays to the given input file.

        Modification history
            20171121-A_schl_ma: written
        """

        # Get variable array
        self.var = self.var_nc[:]

        # Get time array and information
        if (self.time_available):
            self.time = self.time_nc[:]
            self.time_n = len(self.time)
            try:
                self.calendar = self.time_nc.calendar
                self.time_units = self.time_nc.units
            except:
                raise AttributeError("Invalid input: time variable does " + \
                                     "contain calendar or units information")

        # Get lat/lon arrays
        if (self.lat_available):
            self.lat = self.lat_nc[:]
            self.lat_n = len(self.lat)
        if (self.lon_available):
            self.lon = self.lon_nc[:]
            self.lon_n = len(self.lon)

    ###########################################################################

    def check_latlon_available(self):
        """
        Arguments
            None

        Return value
            None

        Description
            Check if latitude and longitude data is available and valid.

        Modification history
            20171121-A_schl_ma: written
        """

        if not (self.lat_available and self.lon_available):
            raise RuntimeError("No lat/lon data for this calculation " + \
                               "available")

    ###########################################################################

    def check_time_available(self):
        """
        Arguments
            None

        Return value
            None

        Description
            Check if time data is available and valid.

        Modification history
            20171121-A_schl_ma: written
        """

        if not (self.time_available):
            raise RuntimeError("No time data for this calculation available")

    ###########################################################################

    def select_region(self, region="global", reset=True):
        """
        Arguments
            region : Region which should be selected
                     Possible values: 'global' or 2D array of the form
                     [[lat_min, lat_max], [lon_min, lon_max]]
            reset  : True: perform calculation on the data of the input
                           netCDF file
                     False: perform calculation on current internal arrays
                            (may be modified by former calculations)

        Return value
            None

        Description
            Shrink down the given arrays to a specified region

        Modification history
            20171121-A_schl_ma: written
        """

        # Check if lat/lon data is available
        self.check_latlon_available()

        # Check if all arguments are valid
        try:
            valid_region_str = True
            valid_region = True
            if (type(region) == str):
                valid_region_str = (region == "global")
            else:
                region = np.array(region)
                valid_region = all([np.shape(region) == (2,2),
                                    region[0,0]>=-90.0, region[0,0]<=90.0,
                                    region[0,1]>=-90.0, region[0,1]<=90.0,
                                    region[1,0]>=0.0, region[1,0]<=360.0,
                                    region[1,1]>=0.0, region[1,1]<=360.0,
                                    region[0,0] < region[0,1],
                                    region[1,0] < region[1,1]])
            valid_reset = (reset is True or reset is False)
        except:
            raise TypeError("Invalid input")
        if (not valid_region_str):
            raise TypeError("Invalid input: region only accepts the " + \
                            "string 'global'")
        if (not valid_region):
            raise ValueError("Invalid input: shape of region array has to " + \
                             "be [[lat_min, lat_max], [lon_min, lon_max]] " + \
                             "with values -90<=lat<=+90, 0<=lon<=360, " + \
                             "lat_min<lat_max and lon_min<lon_max")
        if (not valid_reset):
            raise ValueError("Invalid input: rest has to be True or False")

        # Reset arrays if desired
        if (reset is True):
            self.reset()

        # Get correct region
        if (type(region) == str):
            self.reset()
        else:
            lat_indices = (np.where((self.lat>=region[0,0]) & \
                                    (self.lat<=region[0,1])))[0]
            lon_indices = (np.where((self.lon>=region[1,0]) & \
                                    (self.lon<=region[1,1])))[0]
            if (len(lat_indices) == 0):
                raise ValueError("Invalid input: no data in selected " + \
                                 "region (latitude) available")
            if (len(lon_indices) == 0):
                raise ValueError("Invalid input: no data in selected " + \
                                 "region (longitude) available")
            lat_indices = slice(lat_indices[0], lat_indices[-1]+1)
            lon_indices = slice(lon_indices[0], lon_indices[-1]+1)

            # Try to reshape arrays
            try:
                self.var = self.var[:, lat_indices, lon_indices]
                self.lat = self.lat[lat_indices]
                self.lon = self.lon[lon_indices]
            except:
                raise RuntimeError("Could not select the region from the " + \
                                   "current data, maybe you need to reset " + \
                                   "first")
            self.lat_n = len(self.lat)
            self.lon_n = len(self.lon)

    ###########################################################################

    def gridcell_area(self, lat_bot, lat_top, lon_delta):
        """
        Arguments
            lat_bot   : bottom boundary of the cell [deg]
            lat_top   : top boundary of the cell [deg]
            lon_delta : width of the cell [deg]

        Return value
            The area of the cell [m^2]

        Description
            Calculates the area of a grid cell on the surface of the Earth.

        Modification history
            20171110-A_schl_ma: written
        """

        # Check if all arguments are valid
        try:
            lat_bot = float(lat_bot)
            lat_top = float(lat_top)
            lon_delta = float(lon_delta)
        except:
            raise TypeError("Invalid input")
        valid_input = all([lat_bot>=-90.0, lat_bot<=90.0,
                           lat_top>=-90.0, lat_top<=90.0,
                           lon_delta>=-360.0, lon_delta<=360.0])
        if (not valid_input):
            raise ValueError("Invalid input: Value(s) out of bonds")

        # Get Earth's radius and conversion factor from [deg] to [rad]
        R = Constant.Earth_radius
        conv = np.pi / 180.0

        return abs(R**2 * lon_delta * conv * \
                   (np.sin(lat_top*conv) - np.sin(lat_bot*conv)))

    ###########################################################################

    def map_area(self):
        """
        Arguments
            None

        Return value
            2D numeric array containing the area of each grid cell
            Shape: (latitude, longitude)

        Description
            Calculates the area of each cell of the grid and saves it in the
            member variable

        Modification history
            20171113-A_schl_ma: written
        """

        # Check if lat/lon data is available
        self.check_latlon_available()

        # Calculation is only possible if more grid contains enough cells
        if (self.lat_n<2 or self.lon_n<2):
            raise RuntimeError("Calculation of areas is not possible " + \
                               "because there's not enough data")

        # Calculate latitude interfaces
        lat_interfaces = np.zeros(self.lat_n + 1)

        lat_bottom = (3*self.lat[0] - self.lat[1]) / 2.0
        lat_top = (3*self.lat[self.lat_n-1] - self.lat[self.lat_n-2]) / 2.0

        lat_interfaces[0] = max([lat_bottom, -90.0])
        for i in xrange(1, self.lat_n):
            lat_interfaces[i] = (self.lat[i] + self.lat[i-1]) / 2.0
        lat_interfaces[self.lat_n] = min([lat_top, 90.0])

        # Calculate longitude interfaces
        lon_interfaces = np.zeros(self.lon_n + 1)

        lon_interfaces[0] = (3*self.lon[0] - self.lon[1]) / 2.0
        for i in xrange(1, self.lon_n):
            lon_interfaces[i] = (self.lon[i] + self.lon[i-1]) / 2.0
        lon_interfaces[self.lon_n] = (3*self.lon[self.lon_n-1] -
                                      self.lon[self.lon_n-2]) / 2.0

        # Calculate areas
        map = np.zeros((self.lat_n, self.lon_n))
        for lat_index in xrange(self.lat_n):
            for lon_index in xrange(self.lon_n):
                lon_delta = lon_interfaces[lon_index+1] - \
                            lon_interfaces[lon_index]
                map[lat_index, lon_index] = self.gridcell_area(
                    lat_interfaces[lat_index],
                    lat_interfaces[lat_index + 1],
                    lon_delta)

        return map

    ###########################################################################

    def spatial_average(self, axis="all", weighting=True, region="global",
                        reset=True):
        """
        Arguments
            axis      : Axis along which the average is computed
                        Possible values: 'all', 'lat', 'lon'
            weighting : Weight the variables with the area of the corresponding
                        grid cell
                        Possible values: True or False
            region    : Area for which the average is computed
                        Possible values: 'global' or 2D array of the form
                        [[lat_min, lat_max], [lon_min, lon_max]]
            reset     : True: perform calculation on the data of the input
                              netCDF file
                        False: perform calculation on current internal arrays
                               (may be modified by former calculations)

        Return value
            Averaged value(s) of the given variable

        Description
            Calculates the (weighted) average of a certain variable over a
            specified spatial domain.

        Modification history
            20171114-A_schl_ma: written
        """

        # Check if lat/lon data is available
        self.check_latlon_available()

        # Check if arguments (axis and weighting) are valid
        try:
            valid_axis = any([axis == "lat", axis == "lon", axis == "all"])
            valid_weighting = (weighting is True or weighting is False)
        except:
            raise TypeError("Invalid input: axis or weighting")
        if (not valid_axis):
            raise ValueError("Invalid input: axis has to be 'lat', 'lon' " + \
                             "or 'all'")
        if (not valid_weighting):
            raise ValueError("Invalid input: weighting has to be True or " + \
                             "False")

        # Get correct region (includes check if region and reset are valid)
        self.select_region(region, reset)

        # Get correct weighting
        if (weighting):
            weights = self.map_area()
            print(self.var.shape)
            print(weights.shape)
        else:
            weights = None

        # Average over defined axis for every time step and save member values
        if (axis == "all"):
            var = np.zeros(self.time_n)
            for time_index in xrange(self.time_n):
                var[time_index] = np.average(self.var[time_index],
                                             axis=None, weights=weights)
            self.var = var
            self.lat_available = False
            self.lon_available = False
        elif (axis == "lat"):
            if (weights is not None):
                weights = np.mean(weights, axis=1)
            self.var = np.average(self.var, axis=1, weights=weights)
            self.lat_available = False
        elif (axis == "lon"):
            if (weights is not None):
                weights = np.mean(weights, axis=0)
            self.var = np.average(self.var, axis=2, weights=weights)
            self.lon_available = False

        return self.var

    ###########################################################################

    def temporal_average(self, period="annual", reset=True):
        """
        Arguments
            type  : Defines the period of time over which the average is
                    performed
                    Possible values: 'monthly', 'annual', 'total'
            reset : True: perform calculation on the data of the input
                          netCDF file
                    False: perform calculation on current internal arrays
                           (may be modified by former calculations)

        Return value
            Averaged value(s) of the given variable

        Description
            Calculates the average of a certain variable over a specified
            period of time.

        Modification history
            20171120-A_schl_ma: written
        """

        # Check if time data is available
        self.check_time_available()

        # Check if all arguments are valid
        try:
            valid_period = any([period == "monthly", period == "annual",
                                period == "total"])
        except:
            raise TypeError("Invalid input")
        if (not valid_period):
            raise ValueError("Invalid input: period has to be 'monthly', " + \
                             "'annual' or 'total'")

        # Reset arrays if desired
        if (reset is True):
            self.reset()

        # Average over defined period of time
        if (period == "total"):
            pass
        dates = nc.num2date(self.time, self.time_units, self.calendar)
        months = [d.month for d in dates]
        for m in month:
            pass
        return months

    ###########################################################################
