"""
###############################################################################
# GENERAL ROUTINES FOR OPERATIONS ON A (RECTILINEAR) GRID
###############################################################################
Please consider using of extending existing routines before adding new ones.
Check the header of each routine for documentation.

Contents:
   function gridcell_area

###############################################################################
"""


from constants import Constant

import numpy as np


###############################################################################

def gridcell_area(lon1, lon2, lat1, lat2):
    """
    This function calculates the area of a specified grid cell.

    Arguments
        lon1 : left limit of the cell [deg]
        lon2 : right limit of the cel  [deg]
        lat1 : upper limit of the cell [deg]
        lat2 : lower limit of the cell [deg]

    Return value
        The area of the cell [m^2]

    Description
        Calculates the area of a grid cell on the surface of the Earth.

    Modification history:
        20171110-A_schl_ma: written
    """

    # Check correct ranges of arguments
    try:
        lon1 = float(lon1)
        lon2 = float(lon2)
        lat1 = float(lat1)
        lat2 = float(lat2)
    except:
        raise TypeError("In module 'latlon': in function 'gridcell_area': " + \
                        "Invalid input, function only accepts integers/floats")

    valid_input = (lon1>=0.0 and lon1<=360.0) and \
                  (lon2>=0.0 and lon2<=360.0) and \
                  (lat1>=-90.0 and lat1<=90.0) and \
                  (lat2>=-90.0 and lat2<=90.0)
    if (not valid_input):
        raise ValueError("In module 'latlon': in function gridcell_" + \
                         "area: Latitude/longitude value out of valid bounds")

    # Get Earth's radius and conversion factor from [deg] to [rad]
    R = Constant.Earth_radius
    conv = np.pi / 180.0

    return abs(R**2 * (lon2-lon1) * conv *
               (np.sin(lat2*conv) - np.sin(lat1*conv)))

###############################################################################
