# -*- coding: utf-8 -*-

"""
###############################################################################
GENERAL ROUTINES FOR OPERATIONS ON A (RECTILINEAR) GRID
###############################################################################

Description
    Collection of useful operations on a grid. Please consider extending
    existing routines before adding new ones. Check the header of each routine
    for documentation.

Contents
   function gridcell_area
   function map_area
   function spatial_average

###############################################################################
"""


from constants import Constant

import numpy as np


###############################################################################

def gridcell_area(lat_bot, lat_top, lon_delta):
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

###############################################################################

def map_area(lat, lon):
    """
    Arguments
        lat : 1D monotonic numeric array with at least two elements
              describing the latitudes of the grid
        lon : 1D monotonic numeric array with at least two elements
              describing the longitudes of the grid

    Return value
        2D numeric array containing the area of each grid cell
        Shape: (latitude, longitude)

    Description
        Calculates the area of each cell of a specified grid.

    Modification history
        20171113-A_schl_ma: written
    """

    # Check if all arguments are valid
    try:
        lat = np.array(lat)
        lon = np.array(lon)
        valid_input = all([np.array_equal(lat, np.hstack(lat)),
                           np.array_equal(lon, np.hstack(lon)),
                           len(lat)>1, len(lon)>1])
        monoton_input = (np.array_equal(lat, sorted(lat)) or \
                         np.array_equal(lat, sorted(lat, reverse=True))) and \
                        (np.array_equal(lon, sorted(lon)) or \
                         np.array_equal(lon, sorted(lon, reverse=True)))
    except:
        raise TypeError("Invalid input")

    if (not valid_input):
        raise ValueError("Invalid input: arrays have to be " + \
                         "one-dimensional and contain at least two elements")
    if (not monoton_input):
        raise ValueError("Invalid input: arrays need to be monotonic")

    # Calculate latitude interfaces
    lat_n = len(lat)
    lat_interfaces = np.zeros(lat_n + 1)

    lat_bottom = (3*lat[0] - lat[1]) / 2.0
    lat_top = (3*lat[lat_n-1] - lat[lat_n-2]) / 2.0

    lat_interfaces[0] = max([lat_bottom, -90.0])
    for i in xrange(1, lat_n):
        lat_interfaces[i] = (lat[i] + lat[i-1]) / 2.0
    lat_interfaces[lat_n] = min([lat_top, 90.0])

    # Calculate longitude interfaces
    lon_n = len(lon)
    lon_interfaces = np.zeros(lon_n + 1)

    lon_interfaces[0] = (3*lon[0] - lon[1]) / 2.0
    for i in xrange(1, lon_n):
        lon_interfaces[i] = (lon[i] + lon[i-1]) / 2.0
    lon_interfaces[lon_n] = (3*lon[lon_n-1] - lon[lon_n-2]) / 2.0

    # Calculate areas
    areas = np.zeros((lat_n, lon_n))
    for lat_index in xrange(lat_n):
        for lon_index in xrange(lon_n):
            lon_delta = lon_interfaces[lon_index+1] - lon_interfaces[lon_index]
            areas[lat_index, lon_index] = gridcell_area(
                lat_interfaces[lat_index],
                lat_interfaces[lat_index + 1],
                lon_delta)

    return areas

###############################################################################

def spatial_average(field, lat, lon, axis="all", weights=None):
    """
    Arguments
        field   : 2D numeric array containing the variable which should be
                  averaged
                  Shape: (latitude, longitude)
        lat     : 1D monotonic numeric array with at least two elements
                  describing the latitudes of the grid
        lon     : 1D monotonic numeric array with at least two elements
                  describing the longitudes of the grid
        axis    : Axis along which the average is computed
                  Possible values: 'lat', 'lon', 'all'
        weights : 2D numeric array of weights with the same shape as field
                  'None': Equal weights

    Return value
        Average value(s) of the given variable

    Description
        Calculates the (weighted) average of a certain variable over a
        specified spatial domain.

    Modification history
        20171114-A_schl_ma: written
    """

    # Check if all arguments are valid
    try:
        field = np.array(field)
        lat = np.array(lat)
        lon = np.array(lon)
        valid_field = (np.ndim(field) == 2 and \
                       np.shape(field) == (len(lat), len(lon)))
        valid_latlon = all([np.array_equal(lat, np.hstack(lat)),
                            np.array_equal(lon, np.hstack(lon)),
                            len(lat)>1, len(lon)>1])
        monoton_latlon = (np.array_equal(lat, sorted(lat)) or \
                          np.array_equal(lat, sorted(lat, reverse=True))) and \
                         (np.array_equal(lon, sorted(lon)) or \
                          np.array_equal(lon, sorted(lon, reverse=True)))
        valid_axis = any([axis == "lat", axis == "lon", axis == "all"])
        if (np.array_equal(weights, None)):
            valid_weights = True
        else:
            valid_weights = (np.ndim(weights) == 2 and \
                             np.shape(weights) == (len(lat), len(lon)))
    except:
        raise TypeError("Invalid input")

    if (not valid_field):
        raise ValueError("Invalid input: field array has to be " + \
                         "two-dimensional and match the shape of the " + \
                         "given grid")
    if (not valid_latlon):
        raise ValueError("Invalid input: latitude/longitude arrays have " + \
                         "to be one-dimensional and contain at least two " + \
                         "elements")
    if (not monoton_latlon):
        raise ValueError("Invalid input: latitude/longitude arrays need to " +\
                         "be monotonic")
    if (not valid_axis):
        raise ValueError("Invalid input: axis need to be 'lat', 'lon' or " + \
                         "'all'")
    if (not valid_weights):
        raise ValueError("Invalid input: weights array has to be " + \
                         "two-dimensional and match the shape of the " + \
                         " given grid")

    # Get correct axis index
    axis_index = None
    if (axis == "lat"):
        axis_index = 0
    if (axis == "lon"):
        axis_index = 1

    # Get desired average value
    return np.average(field, axis=axis_index, weights=weights)

###############################################################################
