"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import iris
import os
import c3s_511_utils

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
        
        # config
        self.__project_info__ = dict() # empty project info
        self.__plot_dir__ = '.' + os.sep # default plot directory
        self.__work_dir__ = '.' + os.sep # default work dir

        self.__varname__ = 'var'  # default value
        self.__output_type__ = 'png'  # default ouput file type
        self.__regions__ = {"example":(10,20,-10,-20)} # default regions
        
        self.__verbosity_level__ = 0
        self.__debuginfo__ = "No debug info"
        self.__config__ = dict()
        
        
        # for metadata
        self.__basetags__ = []
        self.__infiles__ = []
        self.authors = ["A_muel_bn", "A_hass_bg", "A_laue_ax", \
                        "A_broe_bj","tbd"] #TODO fill in
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
        
    def __getInfoFromFile__(filename):
        """
        routine to read cfg file
        """
        f = open(filename)
        self.__config__ = imp.load_source('cfg', '', f)
        f.close()
        
