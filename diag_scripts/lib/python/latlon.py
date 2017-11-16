# -*- coding: utf-8 -*-

"""
###############################################################################
GENERAL ROUTINES FOR OPERATIONS ON A (RECTILINEAR) GRID
###############################################################################
"""


from constants import Constant

import numpy as np


class GridOperations(object):
    """
    Description
        Class containing useful operations on a grid. Please consider extending
        existing routines before adding new ones. Check the header of each
        method for documentation.

    Contents
        method gridcell_area
        method map_area
        method spatial_average
    """

    ###########################################################################

    def __init__(self, field, time, lat, lon):
        """
        Arguments
            field : 3D numeric array  describing a certain variable on a grid
                    over a certain time series
                    Shape: (time, latitude, longitude)
            time  : 1D monotonic numeric array describing the temporal progress
                    of the grid
            lat   : 1D monotonic numeric array with at least two elements
                    describing the latitudes of the grid
            lon   : 1D monotonic numeric array with at least two elements
                    describing the longitudes of the grid

        Description
            Initializes class instances.

        Modification history
            20171114-A_schl_ma: written
        """

        # Check if all arguments are valid
        try:
            self.field = np.array(field)
            self.time = np.array(time)
            self.lat = np.array(lat)
            self.lon = np.array(lon)
            self.time_n = len(self.time)
            self.lat_n = len(self.lat)
            self.lon_n = len(self.lon)

            valid_field = (np.ndim(self.field) == 3 and \
                           np.shape(self.field) == \
                                (self.time_n, self.lat_n, self.lon_n))
            valid_time = np.array_equal(self.time, np.hstack(self.time))
            valid_latlon = all([np.array_equal(self.lat, np.hstack(self.lat)),
                                np.array_equal(self.lon, np.hstack(self.lon)),
                                self.lat_n>1, self.lon_n>1])
            mon_time = (np.array_equal(self.time, sorted(self.time)) or \
                          np.array_equal(self.time,
                                         sorted(self.time, reverse=time)))
            mon_latlon = (np.array_equal(self.lat, sorted(self.lat)) or \
                          np.array_equal(self.lat,
                                         sorted(self.lat, reverse=True))) and \
                         (np.array_equal(self.lon, sorted(self.lon)) or \
                          np.array_equal(self.lon,
                                         sorted(self.lon, reverse=True)))
        except:
            raise TypeError("Invalid input")

        if (not valid_field):
            raise ValueError("Invalid input: field array has to be " + \
                             "two-dimensional and match the shape of the " + \
                             "given grid")
        if (not valid_time):
            raise ValueError("Invalid input: time array has to be " + \
                             "one-dimensional")
        if (not valid_latlon):
            raise ValueError("Invalid input: latitude/longitude arrays " + \
                             "have to be one-dimensional and contain at " + \
                             "least two elements")
        if (not mon_time):
            raise ValueError("Invalid input: time array has to be monotonic")
        if (not mon_latlon):
            raise ValueError("Invalid input: latitude/longitude arrays " + \
                             "need to be monotonic")

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

    def spatial_average(self, axis="all", weighting=True):
        """
        Arguments
            axis      : Axis along which the average is computed
                        Possible values: 'all', 'lat', 'lon'
            weighting : Weight the variables with the area of the corresponding
                        grid cell
                        Possible values: True or False

        Return value
            Averaged value(s) of the given variable

        Description
            Calculates the (weighted) average of a certain variable over a
            specified spatial domain.

        Modification history
            20171114-A_schl_ma: written
        """

        # Check if all arguments are valid
        try:
            valid_axis = any([axis == "lat", axis == "lon", axis == "all"])
            valid_weighting = (weighting == True or weighting == False)
        except:
            raise TypeError("Invalid input")

        if (not valid_axis):
            raise ValueError("Invalid input: axis need to be 'lat', 'lon' " + \
                             "or 'all'")
        if (not valid_weighting):
            raise ValueError("Invalid input: weighting has to be True or " + \
                             "False")

        # Get correct weighting
        if (weighting):
            weights = self.map_area()
        else:
            weights = None

        # Average over defined axis for every time step
        if (axis == "all"):
            average = np.zeros(self.time_n)
            for time_index in xrange(self.time_n):
                average[time_index] = np.average(self.field[time_index],
                                                 axis=None, weights=weights)
        elif (axis == "lat"):
            if (weights is not None):
                weights = np.mean(weights, axis=1)
            average = np.average(self.field, axis=1, weights=weights)
        elif (axis == "lon"):
            if (weights is not None):
                weights = np.mean(weights, axis=0)
            average = np.average(self.field, axis=2, weights=weights)

        return average

    ###########################################################################
