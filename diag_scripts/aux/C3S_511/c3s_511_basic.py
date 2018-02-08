"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import iris
import os
import sys

[sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]
import c3s_511_util as utils
from customErrors import *
import warnings
from get_metadata_to_rst import do_report as report
from plot2D import Plot2D


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
        self.__project_info__ = dict()  # empty project info
        self.__plot_dir__ = '.' + os.sep  # default plot directory
        self.__work_dir__ = '.' + os.sep  # default work dir

        self.__varname__ = 'var'  # default value
        self.__output_type__ = 'png'  # default ouput file type
        self.__regions__ = {"example": (10, 20, -10, -20)}  # default regions

        self.__verbosity_level__ = 0  # default information during runtime
        self.__debug_info__ = "No debug info"  # default debug information
        self.__config__ = dict()  # default configuration input

        # for metadata
        self.__basetags__ = []
        self.__infiles__ = []
        self.authors = ["A_muel_bn", "A_hass_bg", "A_laue_ax",
                        "A_broe_bj", "tbd"]  # TODO fill in
        self.diagname = "c3s_511_basic.py"

        self.data = None

    def set_info(self):
        # raise ImplementationError("set_info","This method has to be implemented.")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def read_data(self):
        self.__file_check__()

        return

    def run_diagnostic(self):
        self.__do_overview__()
        self.__do_mean_var__()
        self.__do_trends__()
        self.__do_extremes__()
        self.__do_maturity_matrix__()
        self.__do_gcos_requirements__()

    def __file_check__(self):
        warnings.warn("Implementation Warning", UserWarning)
        return

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
        self.__file_check__()

        datadir = os.getenv('EVT_DATASTORE')
        if datadir is None:
            raise ConfigurationError("Basic_Diagnostic.__init__", "Environment Variable EVT_DATASTORE is not set")
        if os.path.isdir(datadir):
            self.data = iris.load_cube(
                os.path.join(datadir, "test.nc"))
        else:
            raise PathError("Basic_Diagnostic.__init__", "Environment Variable EVT_DATASTORE is not set to valid path.")

        return

    def set_info(self, E, model, var, ref_file, mod_file, cfg):
        """
        gather information for diagnostic
        """

        # TODO: Repair based on E

    #        setags = self._basetags + \
    #            [x.strip() for x in project_info.get('GLOBAL')['tags']]
    #
    #        self._basetags = self._basetags + [self.var]
    #
    #        if var == 'sm':
    #            self._vartype = 'soil moisture'
    #            self._ref_file = ref_file
    #        elif var == 'ts' or var == 'tos':
    #            self._vartype = 'sea surface temperature'
    #            self._ref_file = ref_file
    #        elif var == 'shrubNtreeFrac':
    #            self._vartype = 'shrub and tree'
    #            self._ref_file = ref_file
    #        elif var == 'baresoilFrac':
    #            self._vartype = 'bare soil'
    #            self._self.__config__=utils.__getInfoFromFile__()
    #
    #        self._project_info = project_info
    #
    ## A_laue_ax+
    ##        self._plot_dir=project_info.get('GLOBAL')['plot_dir']
    #        E = ESMValProject(project_info)
    #        plot_dir = E.get_plot_dir()
    #        diag_script = E.get_diag_script_name()
    #        plot_dir = plot_dir + os.sep + diag_script + os.sep
    #        E.ensure_directory(plot_dir)
    #        self._plot_dir = plot_dir
    #        self.E = E
    ## A_laue_ax-
    #
    #        self._work_dir = project_info.get('GLOBAL')['wrk_dir']
    #        self._climo_dir = project_info.get('GLOBAL')['climo_dir']
    #        self._mod_line = model.split_entries()
    #        self._field_type = project_info['RUNTIME']['derived_field_type']
    #        self.var = var
    #                   ['currDiag'].variables]
    #        whichvar = [i for i, x in enumerate(allvars) if x][0]
    #        self._ref_type = project_info['RUNTIME']['currDiag']
    #        self._ref_type = self._ref_type.variables[whichvar].ref_model
    #        self._mod_type = project_info['RUNTIME']['project']
    #
    #        self.output_type = project_info['GLOBAL']['output_file_type']
    #
    #        self._baref_file = ref_file
    #
    #            if self.cfg.cmap_globmeants[-3:-1] == "_r":
    #                self.cfg.cmap_globmeants = self.cfg.cmap_globmeants[:-2]
    #            else:
    #                self.cfg.cmap_globmeants = self.cfg.cmap_globmeants + "_r"
    #
    #        elif var == 'grassNcropFrac':
    #            self._vartype = 'grass and crop'
    #            self._ref_file = ref_file
    #        elif var == 'alb':
    #            self._vartype = 'albedo'
    #            self._ref_file = ref_file
    #        elif var == 'fAPAR':
    #            self._vartype = 'fAPAR'
    #            self._ref_file = ref_file
    #        elif var == 'burntArea':
    #            self._vartype = 'burned area'
    #            self._ref_file = ref_file
    #        else:
    #            assert False, 'This variable is not implemented yet!'
    #
    #        self._mod_file = mod_file
    #
    #        self._get_output_rootname()
    #
    #        if self.cfg.regionalization:
    #            self._reg_file = "./diag_scripts/aux" + \
    #                "/LMU_ESACCI-diagnostics/Shapefiles" + \
    #                os.sep + self.cfg.shape
    #            self._load_regionalization_shape()

    def __do_mean_var__(self):

        Plot2D(self.slice)

        return

    def write_reports(self):

        print report
        report(["/media/bmueller/Work/ESMVAL_res/work/reports/sphinx/source/latlon.png"], "test")
        return
