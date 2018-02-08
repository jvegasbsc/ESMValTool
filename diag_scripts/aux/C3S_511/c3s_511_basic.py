"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import iris
import os
import sys
import numpy as np

[sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]
import c3s_511_util as utils
from customErrors import ConfigurationError, PathError, EmptyContentError
import warnings
from get_metadata_to_rst import do_report as report
from plot2D import Plot2D
from esmval_lib import ESMValProject

import inspect as INS


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

class __Diagnostic_skeleton__(object):
    """
    Basic class to implement any kind of diagnostic
    """

    def __init__(self, **kwargs):
        super(__Diagnostic_skeleton__, self).__init__(**kwargs)
        """
        Default values to experiment with the diagnostics
        """

        # config
        self.__project_info__ = None  # empty project info
        self.__plot_dir__ = '.' + os.sep  # default plot directory
        self.__work_dir__ = '.' + os.sep  # default work dir

        self.__varname__ = 'var'  # default value
        self.__output_type__ = 'png'  # default ouput file type
        self.__regions__ = {"example": (10, 20, -10, -20)}  # default regions

        self.__verbosity_level__ = 0  # default information during runtime
        self.__debug_info__ = "No debug info"  # default debug information
        self.__config__ = dict()  # default configuration input

        self.__basetags__ = []
        self.__infile__ = None
        self.__inpath__= None
        self.report_dict = dict()
        
        self.authors = ["A_muel_bn", "A_hass_bg", "A_laue_ax",
                        "A_broe_bj", "tbd"]  # TODO fill in
        self.diagname = "c3s_511_basic.py"

        self.sp_data = None

    def set_info(self, **kwargs):
        # raise ImplementationError("set_info","This method has to be implemented.")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def read_data(self):
        warnings.warn("Implementation Warning", UserWarning)

        return

    def run_diagnostic(self):
        self.__do_overview__()
        self.__do_mean_var__()
        self.__do_trends__()
        self.__do_extremes__()
        self.__do_maturity_matrix__()
        self.__do_gcos_requirements__()

    def __do_overview__(self):
        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)

        return

    def __do_mean_var__(self):
        
        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)

        return

    def __do_trends__(self):
        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)

        return

    def __do_extremes__(self):
        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)

        return

    def __do_maturity_matrix__(self):
        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)

        return

    def __do_gcos_requirements__(self):
        self.__prepare_report__()

        warnings.warn("Implementation Warning", UserWarning)
        return

    def __prepare_report__(self):
        warnings.warn("Implementation Warning", UserWarning)
        return

    def write_reports(self):
        warnings.warn("Implementation Warning", UserWarning)
        return


class Basic_Diagnostic(__Diagnostic_skeleton__):
    """    def write_reports(self):

        return
    class to implement basic diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def __init__(self, **kwargs):
        super(Basic_Diagnostic, self).__init__(**kwargs)

        #        self.__config__ = utils.__getInfoFromFile__("")

    #        self.data = iris.load_cube("/media/bmueller/Work/ESMVAL_res/work/climo/CMIP5/CMIP5_Amon_historical_MPI-ESM-P_r1i1p1_T2Ms_ts_1991-2005.nc")
    #        self.slice = self.data.collapsed("time",iris.analysis.MEAN)

    # self.slice = iris.load_cube('/media/bmueller/Work/GIT/ESMValTool-private_base/diag_scripts/aux/C3S_511/plots/test_latlon.nc')

    def read_data(self):

        if os.path.isfile(self.__infile__):
            self.sp_data = iris.load_cube(self.__infile__)
        else:
            self.read_data_mock()
            
            print("self.__infile__ was not found! Generic data used instead.")
#            raise PathError("Basic_Diagnostic.__init__", "self.__infile__ is not set to valid path.")
            
        return

    def read_data_mock(self):
        """
        reads artificial data for testing uses
        TODO: Remove after loading is tested.
        """
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/testdata")
        self.sp_data = iris.load_cube(
            os.path.join(datadir, "test.nc"))

    def set_info(self, **kwargs):
        """
        gather information for diagnostic
        """
        self.__project_info__ = kwargs.get('proj_info', None)
        if not isinstance(self.__project_info__, ESMValProject):
            raise EmptyContentError("proj_info", "Element is empty.")

        self.__plot_dir__ = self.__project_info__.get_plot_dir()
        self.__work_dir__ = self.__project_info__.get_work_dir()
        
        self.__varname__ = self.__project_info__.get_currVars()  # default value
        if not len(self.__varname__)==1:
            raise EmptyContentError("self.__varname__", "Element is of wrong length.")
        self.__output_type__ = 'png'  # default ouput file type
        self.__regions__ = {"example": (10, 20, -10, -20)}  # default regions

        self.__verbosity_level__ = self.__project_info__.get_verbosity() # default information during runtime
        self.__debug_info__ = "No debug info"  # default debug information
        self.__config__ = self.__project_info__.get_configfile()  # default configuration input

        # for metadata
        self.__basetags__ = self.__project_info__.get_all_tags()
        file_info = self.__project_info__.get_all_clim_models(variables=self.__varname__[0])
        self.__infile__ = file_info.keys()
        if not len(self.__infile__)==1:
            raise EmptyContentError("self.__infile__", "Element is of wrong length.")
        else:
            self.__infile__ = self.__infile__[0]
            
        self.__inpath__ = file_info[self.__infile__]["dir"]
        
        
    def __do_overview__(self):
        
        missing_values = self.sp_data.copy()
        missing_values.data = np.ones_like(missing_values.data)
        
        missing_values = missing_values.collapsed("time",iris.analysis.SUM)
        missing_values.data = missing_values.data / \
            float(len(self.sp_data.coord("time").points))
            
        Plot2D(missing_values)
        
        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)

        return

    def __do_mean_var__(self):

        # Plot2D(self.slice)
        return

    def write_reports(self):
        #        print report
        #        report(["/media/bmueller/Work/ESMVAL_res/work/reports/sphinx/source/latlon.png"],"test")
        return
