"""
###############################################################################
# GENERAL ROUTINES FOR OPERATIONS ON A (RECTILINEAR) GRID
###############################################################################

Description:
    Collection of useful operations on a grid. Please consider extending
    existing routines before adding new ones. Check the header of each routine
    for documentation.

Contents:
   function gridcell_area
   function map_area

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

    Modification history:
        20171110-A_schl_ma: written
    """

    # Check correct ranges of arguments
    try:
        lat_bot = float(lat_bot)
        lat_top = float(lat_top)
        lon_delta = float(lon_delta)
    except:
        raise TypeError("Invalid input: No numerical values given")

    valid_input = (lat_bot>=-90.0 and lat_bot<=90.0) and \
                  (lat_top>=-90.0 and lat_top<=90.0) and \
                  (lon_delta>=-360.0 and lon_delta<=360.0)
    if (not valid_input):
        raise ValueError("Invalid input: Value(s) out of bonds")

    # Get Earth's radius and conversion factor from [deg] to [rad]
    R = Constant.Earth_radius
    conv = np.pi / 180.0

    return abs(R**2 * lon_delta * conv * \
               (np.sin(lat_top*conv) - np.sin(lat_bot*conv)))

###############################################################################

def map_area(latitudes, longitudes):
    """
    Arguments
        latitudes  : monotonic numeric array of dimension 1 with at least two
                     elements containing the latitudes of the grid
        longitudes : monotonic numeric array of dimension 1 with at least two
                     elements containing the longitudes of the grid

    Return value
        Numeric array of dimension 2 containing the area of each grid cell
        Shape: (latitude, longitude)

    Description
        Calculates the area of  each cell of a specified grid.

    Modification history:
        20171113-A_schl_ma: written
    """

    # Check correct dimension of the arguments
    valid_input = True
    monotonic_input = True
    try:
        valid_input = (np.ndim(latitudes)==1 and np.ndim(longitudes)==1) and \
                      (len(latitudes)>1 and len(longitudes)>1)
        monotonic_input = (latitudes == sorted(latitudes) or \
                           latitudes == sorted(latitudes, reverse=True)) and \
                          (longitudes == sorted(longitudes) or \
                           longitudes == sorted(longitudes, reverse=True))
    except:
        raise TypeError("Invalid input: function only accepts numeric arrays")

    if (not valid_input):
        raise ValueError("Invalid input: arrays need to be of dimension 1 " + \
                         "and contain at least two elements")
    if (not monotonic_input):
        raise ValueError("Invalid input: arrays need to be monotonic")

    # Calculate latitude interfaces
    lat_n = len(latitudes)
    lat_interfaces = np.zeros(lat_n + 1)

    lat_bottom = (3*latitudes[0] - latitudes[1]) / 2.0
    lat_top = (3*latitudes[lat_n-1] - latitudes[lat_n-2]) / 2.0

    lat_interfaces[0] = max([lat_bottom, -90.0])
    for i in xrange(1, lat_n):
        lat_interfaces[i] = (latitudes[i] + latitudes[i-1]) / 2.0
    lat_interfaces[lat_n] = min([lat_top, 90.0])

    # Calculate longitude interfaces
    lon_n = len(longitudes)
    lon_interfaces = np.zeros(lon_n + 1)

    lon_interfaces[0] = (3*longitudes[0] - longitudes[1]) / 2.0
    for i in xrange(1, lon_n):
        lon_interfaces[i] = (longitudes[i] + longitudes[i-1]) / 2.0
    lon_interfaces[lon_n] = (3*longitudes[lon_n-1] - longitudes[lon_n-2]) / 2.0

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
