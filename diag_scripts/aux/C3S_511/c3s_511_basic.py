"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import iris
import os
import lib.c3s_511_utils as utils

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
    
    
    def set_info(self):
        
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
    
    
    def __prepare_report__(self):self.
        
        return



class Basic_Diagnostic(Diagnostic_skeleton):
    """
    class to implement basic diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def __init__(self, **kwargs):
        super(Basic_Diagnostic, self).__init__(**kwargs)
        
        self.__config__=utils.__getInfoFromFile__()
        
    def set_info(self, project_info, model, var, ref_file, mod_file, cfg):
        """
        gather information for diagnostic
        """

        self.__config__=utils.__getInfoFromFile__()

        self._project_info = project_info

# A_laue_ax+
#        self._plot_dir=project_info.get('GLOBAL')['plot_dir']
        E = ESMValProject(project_info)
        plot_dir = E.get_plot_dir()
        diag_script = E.get_diag_script_name()
        plot_dir = plot_dir + os.sep + diag_script + os.sep
        E.ensure_directory(plot_dir)
        self._plot_dir = plot_dir
        self.E = E
# A_laue_ax-

        self._work_dir = project_info.get('GLOBAL')['wrk_dir']
        self._climo_dir = project_info.get('GLOBAL')['climo_dir']
        self._mod_line = model.split_entries()
        self._field_type = project_info['RUNTIME']['derived_field_type']
        self.var = var
        allvars = [x.var == var for x in project_info['RUNTIME']
                   ['currDiag'].variables]
        whichvar = [i for i, x in enumerate(allvars) if x][0]
        self._ref_type = project_info['RUNTIME']['currDiag']
        self._ref_type = self._ref_type.variables[whichvar].ref_model
        self._mod_type = project_info['RUNTIME']['project']

        self.output_type = project_info['GLOBAL']['output_file_type']

        self._basetags = self._basetags + \
            [x.strip() for x in project_info.get('GLOBAL')['tags']]

        self._basetags = self._basetags + [self.var]

        if var == 'sm':
            self._vartype = 'soil moisture'
            self._ref_file = ref_file
        elif var == 'ts' or var == 'tos':
            self._vartype = 'sea surface temperature'
            self._ref_file = ref_file
        elif var == 'shrubNtreeFrac':
            self._vartype = 'shrub and tree'
            self._ref_file = ref_file
        elif var == 'baresoilFrac':
            self._vartype = 'bare soil'
            self._ref_file = ref_file

            if self.cfg.cmap_globmeants[-3:-1] == "_r":
                self.cfg.cmap_globmeants = self.cfg.cmap_globmeants[:-2]
            else:
                self.cfg.cmap_globmeants = self.cfg.cmap_globmeants + "_r"

        elif var == 'grassNcropFrac':
            self._vartype = 'grass and crop'
            self._ref_file = ref_file
        elif var == 'alb':
            self._vartype = 'albedo'
            self._ref_file = ref_file
        elif var == 'fAPAR':
            self._vartype = 'fAPAR'
            self._ref_file = ref_file
        elif var == 'burntArea':
            self._vartype = 'burned area'
            self._ref_file = ref_file
        else:
            assert False, 'This variable is not implemented yet!'

        self._mod_file = mod_file

        self._get_output_rootname()

        if self.cfg.regionalization:
            self._reg_file = "./diag_scripts/aux" + \
                "/LMU_ESACCI-diagnostics/Shapefiles" + \
                os.sep + self.cfg.shape
            self._load_regionalization_shape()
