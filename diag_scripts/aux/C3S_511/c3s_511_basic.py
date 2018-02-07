"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import numpy as np
import os
import pdb
# import matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#import matplotlib.dates as mdates
# from netCDF4 import Dataset

# global installation
from geoval.core.data import GeoData
from geoval.core.mapping import SingleMap
import extended_data
import extended_mapping
from easyhov import easyhov_diff
from esmval_lib import ESMValProject
from ESMValMD import ESMValMD
# from GeoData_mapping import *

# import ConfigParser
import csv
import imp
import shapefile as shp
import glob
import math
import tempfile
import datetime
# from dateutil.relativedelta import relativedelta
# import subprocess
# import fnmatch

from scipy import stats
from cdo import Cdo

# All packages checked

# ignored GLOBAL values:
# * verbosity_level
# * exit_on_warning
# * debuginfo

# * max_data_filesize
# * max_data_blocksize

# * write_plots
# * write_netcdf
# * write_plot_vars

# * force_processing


class Diagnostic_skeleton(object):
    """
    Basic class to implement any kind of diagnostic
    """

    def __init__(self, **kwargs):
        super(Diagnostic_skeleton, self).__init__(**kwargs)
        """
        Default values to experiment with the diagnostics
        """
#        self._project_info = {}
#        self._mod_type = 'model'
#        self._ref_type = 'reference'
#        self._plot_dir = '.' + os.sep
#        self._work_dir = '.' + os.sep
#
#        self._vartype = 'some variable'  # default value as there must be one
#        self.output_type = 'png'  # default ouput file type
#        self.overview = False
#        self._regions = None
#        self._changed = False
        
        

        self.__basetags = []
        self.__infiles = []

        # additional meta data
        self.authors = ["A_muel_bn", "A_hass_br", "tbd"]
        self.diagname = "c3s_511_basic.py"


    def read_data(self):
        
        self.__file_check__()
        
        return
    
    
    def run_diagnostic(self):
        
        self.__do_overview__()
        self.__do_gcos__()
        self.__do_mean_var__()
        self.__do_trends__()
        
    
    def __file_check__(self):
    
        return
    
    
    def __do_overview__(self):
        
        self.__prepare_report__()
        
        return
    
    
    def __do_mean_var__(self):
        
        self.__prepare_report__()
        
        return
    
    
    def __do_trends__(self):
        
        self.__prepare_report__()
        
        return
    
    
    def __do_extremes__(self):
        
        self.__prepare_report__()
        
        return
    
    
    def __do_maturity_matrix__(self):
        
        self.__prepare_report__()
        
        return
    
    
    def __prepare_report__(self):
        
        return



class Basic_Diagnostic(Diagnostic_skeleton):
    """
    class to implement basic diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def __init__(self, **kwargs):
        super(Basic_Diagnostic, self).__init__(**kwargs)
        
