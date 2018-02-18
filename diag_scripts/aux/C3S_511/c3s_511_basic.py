"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import iris
import os
import sys
import shutil
import errno
import numpy as np
import random, string
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

[sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]
import c3s_511_util as utils
from customErrors import ImplementationError, ConfigurationError, PathError, EmptyContentError
import warnings
from get_metadata_to_rst import do_report as report
from get_metadata_to_rst import do_smm_table
from get_metadata_to_rst import do_gcos_table
from plot import Plot2D, PlotHist
from esmval_lib import ESMValProject
from ESMValMD import ESMValMD


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
                        "A_broe_bj", "A_mass_fr", "A_nico_nd",
                        "A_schl_mn", "A_bock_ls"]  # TODO fill in
        self.diagname = "c3s_511_skeleton.py"

        self.sp_data = None

    def set_info(self, **kwargs):
        raise ImplementationError("set_info","This method has to be implemented.")
        return

    def read_data(self):
        raise ImplementationError("read_data","This method has to be implemented.")
        return

    def run_diagnostic(self):
        self.__do_overview__()
        self.__do_mean_var__()
        self.__do_trends__()
        self.__do_extremes__()
        self.__do_maturity_matrix__()
        self.__do_gcos_requirements__()

    def __do_overview__(self):
#        self.__prepare_report__()
        raise ImplementationError("__do_overview__","This method has to be implemented.")
        return

    def __do_mean_var__(self):
#        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_trends__(self):
#        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_extremes__(self):
#        self.__prepare_report__()
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_maturity_matrix__(self):
#        self.__prepare_report__()
        raise ImplementationError("__do_maturity_matrix__","This method has to be implemented.")
        return

    def __do_gcos_requirements__(self):
#        self.__prepare_report__()
        raise ImplementationError("__do_gcos_requirements__","This method has to be implemented.")
        return

    def __prepare_report__(self):
        warnings.warn("Implementation Warning", UserWarning)
        return

#    def write_reports(self):
#        #currently done with preparation
#        warnings.warn("Implementation Warning", UserWarning)
#        return


class Basic_Diagnostic(__Diagnostic_skeleton__):
    """    def write_reports(self):

        return
    class to implement basic diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def __init__(self, **kwargs):
        super(Basic_Diagnostic, self).__init__(**kwargs)
        self.diagname = "c3s_511_basic.py"


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
        else:
            self.__varname__ = self.__varname__[0]
        self.__output_type__ = 'png'  # default ouput file type for the basic diagnostics
        self.__regions__ = {"example": (10, 20, -10, -20)}  # default regions

        self.__verbosity_level__ = self.__project_info__.get_verbosity() # default information during runtime
        self.__debug_info__ = "No debug info"  # default debug information
        self.__config__ = utils.__getInfoFromFile__(self.__project_info__.get_configfile())  # set configuration input
        
        # for metadata
        self.__basetags__ = self.__project_info__.get_all_tags() + [self.__varname__]
        file_info = self.__project_info__.get_all_clim_models(variables=self.__varname__)
        self.__infile__ = file_info.keys()
        if not len(self.__infile__)==1:
            raise EmptyContentError("self.__infile__", "Element is of wrong length.")
        else:
            self.__infile__ = self.__infile__[0]
            
        self.__inpath__ = file_info[self.__infile__]["dir"]
        self.__time_period__ = "-".join([file_info[self.__infile__]["start_year"],file_info[self.__infile__]["end_year"]])
        self.__dataset_id__ = [file_info[self.__infile__]["name"], file_info[self.__infile__]["case_name"], file_info[self.__infile__]["ensemble"], file_info[self.__infile__]["var"]]
        
        self.__basic_filename__ = "_".join(self.__dataset_id__ + [self.__time_period__])


    def read_data(self):

        if os.path.isfile(self.__infile__):
            self.sp_data = iris.load_cube(self.__infile__)
        else:
            self.__read_data_mock__()
            
            print("self.__infile__ was not found! Generic data used instead.")
#            raise PathError("Basic_Diagnostic.__init__", "self.__infile__ is not set to valid path.")
            
        return


    def __read_data_mock__(self):
        """
        reads artificial data for testing uses
        TODO: Remove after loading is tested.
        """
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tests/testdata")
        self.sp_data = iris.load_cube(
            os.path.join(datadir, "test.nc"))       
        
        
    def __do_overview__(self):
        
        this_file = "overview"
        
        # TODO all the other plots
        
        list_of_plots=[]
        
        # Lat lon plot of fractional available measurements
        missing_values_2D = self.sp_data.copy()
        missing_values_2D.data = np.ma.masked_array(missing_values_2D.data, mask=np.isnan(missing_values_2D.data))
        missing_values_2D.data = missing_values_2D.data*0.+1.
        missing_values_2D = missing_values_2D.collapsed("time", iris.analysis.SUM)
        missing_values_2D.data = missing_values_2D.data / \
            float(len(self.sp_data.coord("time").points))
         
        # plotting routine
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_frac_avail_lat_lon" + "." + self.__output_type__
        list_of_plots.append(filename)
        x=Plot2D(missing_values_2D)
        x.plot(title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")").savefig(filename)
        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_overview'],
                 str('Overview on latitude/longitude availablility of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 '#C3S' + 'fravlalo' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # Lat time plot of fractional available measurements
        missing_values_Lati = self.sp_data.copy()
        missing_values_Lati.data = np.ma.masked_array(missing_values_Lati.data, mask=np.isnan(missing_values_Lati.data))
        missing_values_Lati.data = missing_values_Lati.data*0.+1.
        missing_values_Lati = missing_values_Lati.collapsed("longitude", iris.analysis.SUM)
        missing_values_Lati.data = missing_values_Lati.data / \
            float(len(self.sp_data.coord("longitude").points))
         
        # plotting routine
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_frac_avail_lat_time" + "." + self.__output_type__
        list_of_plots.append(filename)
        x = Plot2D(missing_values_Lati)
        fig = x.plot(summary_plot=True, title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
        fig.savefig(filename)
        plt.close(fig)
        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_overview'],
                 str('Overview on latitude/time availablility of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 '#C3S' + 'fravlati' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # Lon time plot of fractional available measurements
        missing_values_Loti = self.sp_data.copy()
        missing_values_Loti.data = np.ma.masked_array(missing_values_Loti.data, mask=np.isnan(missing_values_Loti.data))
        missing_values_Loti.data = missing_values_Loti.data*0.+1.
        missing_values_Loti = missing_values_Loti.collapsed("latitude", iris.analysis.SUM)
        missing_values_Loti.data = missing_values_Loti.data / \
            float(len(self.sp_data.coord("latitude").points))
         
        # plotting routine
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_frac_avail_lon_time" + "." + self.__output_type__
        list_of_plots.append(filename)
        x=Plot2D(missing_values_Loti)
        fig = x.plot(summary_plot=True, title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
        fig.savefig(filename)
        plt.close(fig)
        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_overview'],
                 str('Overview on longitude/time availablility of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 '#C3S' + 'fravloti' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # Lon time plot of fractional available measurements
        all_data = self.sp_data.copy()
        all_data.data = np.ma.masked_array(all_data.data, mask=np.isnan(all_data.data))
         
        # plotting routine
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_hist_all_vals" + "." + self.__output_type__
        list_of_plots.append(filename)
        x=PlotHist(all_data)
        fig = x.plot(title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
        fig.savefig(filename)
        plt.close(fig)
        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_overview'],
                 str('Full temporal and spatial histogram of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                 '#C3S' + 'histall' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)


        # dimension information
        lon_range = (self.sp_data.coord("longitude").points)
        lat_range = (self.sp_data.coord("latitude").points)
        tim_range = (self.sp_data.coord("time").points)
        
        lon_range_spec = utils.__minmeanmax__(lon_range)
        lat_range_spec = utils.__minmeanmax__(lat_range)
        tim_range_spec = utils.__minmeanmax__(tim_range)


        lon_freq = np.diff(lon_range)
        lat_freq = np.diff(lat_range)
        tim_freq = np.diff(tim_range)
        
        lon_freq_spec = utils.__minmeanmax__(lon_freq)
        lat_freq_spec = utils.__minmeanmax__(lat_freq)
        tim_freq_spec = utils.__minmeanmax__(tim_freq)
        
        overview_dict=collections.OrderedDict()
        
        overview_dict.update({'longitude range [degrees]': collections.OrderedDict([("min",str(lon_range_spec[0])), ("max",str(lon_range_spec[2]))])})
        overview_dict.update({'longitude frequency [degrees]': collections.OrderedDict([("min",str(lon_freq_spec[0])), ("average",str(lon_freq_spec[1])), ("max",str(lon_freq_spec[2]))])})
        overview_dict.update({'latitude range [degrees]': collections.OrderedDict([("min",str(lat_range_spec[0])), ("max",str(lat_range_spec[2]))])})
        overview_dict.update({'latitude frequency [degrees]': collections.OrderedDict([("min",str(lat_freq_spec[0])), ("average",str(lat_freq_spec[1])), ("max",str(lat_freq_spec[2]))])})
        overview_dict.update({'temporal range [days since 01/01/1950]': collections.OrderedDict([("min",str(tim_range_spec[0])), ("max",str(tim_range_spec[2]))])})
        overview_dict.update({'temporal frequency [days]': collections.OrderedDict([("min",str(tim_freq_spec[0])), ("average",str(tim_freq_spec[1])), ("max",str(tim_freq_spec[2]))])})
        
        # produce reports
        self.__prepare_report__(content={"text":overview_dict,"plots":list_of_plots}, filename=this_file.upper())

        return
    
    def __prepare_report__(self,**kwargs):
        
        # TODO specify sphinx structure
        
        content = kwargs.get('content', [])
        if not isinstance(content, (list, dict)):
            raise TypeError("content", "Element is not a list, nor a dict.")
            
        rand_str = lambda n: ''.join([random.choice(string.lowercase) for i in xrange(n)])
            
        filename = kwargs.get('filename', rand_str(10))
        if not isinstance(filename, str):
            raise TypeError("filename", "Element is not a string.")
            
        report(content,filename,self.__work_dir__)
        return
    
    
    def __do_maturity_matrix__(self):
        
        this_file = "System maturity matrix"
        
        expected_input = self.__work_dir__ + os.sep + "smm_input" + os.sep + self.__basic_filename__ + "_smm_expert.csv"
        if not os.path.isfile(expected_input):
            try:
                os.makedirs(os.path.dirname(expected_input))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
            shutil.copy2(os.path.dirname(os.path.realpath(__file__)) + "/lib/predef/empty_smm_expert.csv", expected_input)
            print("*******************************************************************************WARNING*******************************************************************************")
            print("Expected " + this_file + " input file " + expected_input + " not found!")
            print("Created dummy file instead. Please fill in appropriate values and rerun!")
            print("(This won't fail if you do not do so and produce a white matrix!!!!)")
            print("*******************************************************************************WARNING*******************************************************************************")
            captionerror = True
        else:
            print("Processing " + this_file + " input file: " + expected_input) 
            captionerror = False

        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_file.split()) + "." + self.__output_type__
        fig = do_smm_table(expected_input, os.path.dirname(os.path.realpath(__file__)) + "/lib/predef/smm_definitions.csv")
        fig.savefig(filename)
        plt.close(fig)
        
        if captionerror:
            caption = "The expected input file: " + expected_input + " was not found and an empty dummy file created, therefore this plot is blank. Please edit the expected file!"
        else:
            caption = str(this_file + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')

        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_SMM'],
                 caption,
                 '#C3S' + 'SMM' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        self.__prepare_report__(content={"plots":[filename]}, filename="".join(this_file.upper().split()))
        
        return
    
    
    def __do_gcos_requirements__(self):
        
        this_file = "GCOS requirements"
        
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_file.split()) + "." + self.__output_type__
        fig = do_gcos_table(os.path.dirname(os.path.realpath(__file__)) + "/example_csvs/example_gcos_expert.csv", os.path.dirname(os.path.realpath(__file__)) + "/example_csvs/example_gcos_reference.csv")
        fig.savefig(filename)
        plt.close(fig)
        
        caption = str(this_file + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')
        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_GCOS'],
                 caption,
                 '#C3S' + 'GCOS' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        self.__prepare_report__(content={"plots":[filename]}, filename="".join(this_file.upper().split()))
        
        return
