"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import iris
import iris.coord_categorisation
from iris.analysis import stats as cubestats
import iris.coords as icoords
import os
import sys
import shutil
import errno
import numpy as np
import random, string
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import stats
import datetime

[sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]

import c3s_511_util as utils
from customErrors import ImplementationError, ConfigurationError, PathError, EmptyContentError
import warnings
from get_metadata_to_rst import do_report as report
from get_metadata_to_rst import do_smm_table
from get_metadata_to_rst import do_gcos_table
from plot import Plot2D_2, PlotHist, Plot2D_blank
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
        self.__gcos_dict__ = dict()
        # save GCOS requirements
        self.__gcos_dict__.update({"Frequency":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Resolution":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Accuracy":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Stability":{"value":None, "unit":None}})
        
        self.authors = ["A_muel_bn", "A_hass_bg", "A_laue_ax",
                        "A_broe_bj", "A_mass_fr", "A_nico_nd",
                        "A_schl_mn", "A_bock_ls"]  # TODO fill in
        self.diagname = "Diagnostic_skeleton.py"
        self.CDS_ID = "XXX_XX_01"

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
#        self.__do_trends__()
#        self.__do_extremes__()
#        self.__do_sectors__()
#        self.__do_maturity_matrix__()
#        self.__do_gcos_requirements__()
#        self.__do_esmvalidation__()

    def __do_overview__(self):
        self.__do_report__(content={},filename="do_overview_default")
        raise ImplementationError("__do_overview__","This method has to be implemented.")
        return

    def __do_mean_var__(self):
        self.__do_report__(content={},filename="do_mean_var_default")
        raise ImplementationError("__do_mean_var__","This method has to be implemented.")
        return

    def __do_trends__(self):
        self.__do_report__(content={},filename="do_trends_default")
        raise ImplementationError("__do_trends__","This method has to be implemented.")
        return

    def __do_extremes__(self):
        self.__do_report__(content={},filename="do_extremes_default")
        warnings.warn("Implementation Warning", UserWarning)
        return
    
    def __do_sectors__(self):
        self.__do_report__(content={},filename="do_sectors_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_maturity_matrix__(self):
        self.__do_report__(content={},filename="do_maturity_matrix_default")
        raise ImplementationError("__do_maturity_matrix__","This method has to be implemented.")
        return

    def __do_gcos_requirements__(self):
        self.__do_report__(content={},filename="do_gcos_requirements_default")
        raise ImplementationError("__do_gcos_requirements__","This method has to be implemented.")
        return
    
    def __do_esmvalidation__(self):
        self.__do_report__(content={},filename="do_esmvalidation_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_report__(self):
        raise ImplementationError("__do_report__","This method has to be implemented.")
        return

#    def write_reports(self):
#        #currently done with preparation
#        warnings.warn("Implementation Warning", UserWarning)
#        return


class Basic_Diagnostic(__Diagnostic_skeleton__):
    """    
    class to implement basic diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def __init__(self, **kwargs):
        super(Basic_Diagnostic, self).__init__(**kwargs)
        self.diagname = "Basic_Diagnostic.py"


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
        try:
            self.__dataset_id__ = [file_info[self.__infile__]["name"], file_info[self.__infile__]["case_name"], file_info[self.__infile__]["ensemble"], file_info[self.__infile__]["var"]]
        except:
            self.__dataset_id__ = [file_info[self.__infile__]["name"], file_info[self.__infile__]["mip"], file_info[self.__infile__]["experiment"], file_info[self.__infile__]["ensemble"], file_info[self.__infile__]["var"]]

        self.__basic_filename__ = "_".join(self.__dataset_id__ + [self.__time_period__])
        
        self.dimensions = np.array(["time", "latitude", "longitude"]) # TODO: get from cube
        

    def read_data(self):

        if os.path.isfile(self.__infile__):
            self.sp_data = iris.load_cube(self.__infile__)
            self.sp_data.data = np.ma.masked_array(self.sp_data.data, mask=np.isnan(self.sp_data.data))
            self.sp_data.coord('latitude').guess_bounds()
#            if not self.sp_data.coords('grid_longitude'):
#            mima_lon=np.nanpercentile(self.sp_data.coord('longitude').points,[0,100])
#            if np.all(mima_lon <= 360) and np.all(mima_lon>=0):
#                print self.sp_data.coord('longitude')
#                long_2 = self.sp_data.coord('longitude').points.copy()
#                long_2[long_2>180]=long_2[long_2>180]-360
#                long_2_coord = icoords.Coord(long_2, units='degrees', long_name=u'longitude_2', var_name='lon_2')
#                print long_2_coord
#                self.sp_data.add_aux_coord(long_2_coord)
#                print self.sp_data
#                self.sp_data.coord('longitude_2').guess_bounds()
            self.sp_data.coord('longitude').guess_bounds()
            if self.sp_data.units == "no-unit":
                self.sp_data.units = '1'
        else:
            self.__read_data_mock__()
            
            print("self.__infile__ was not found! Generic data used instead.")
#            raise PathError("Basic_Diagnostic.__init__", "self.__infile__ is not set to valid path.")
            
            
        # Human readable time
        t_info = str(self.sp_data.coord("time").units)
        origin = datetime.datetime.strptime("_".join(t_info.split(" ")[-2:]), '%Y-%m-%d_%H:%M:%S')
        if t_info.split(" ")[0] == 'days':
            tim_read = [origin + datetime.timedelta(days=dt) for dt in self.sp_data.coord("time").points]
        elif t_info.split(" ")[0] == 'seconds':
            tim_read = [origin + datetime.timedelta(seconds=dt) for dt in self.sp_data.coord("time").points]
        elif t_info.split(" ")[0] == 'microseconds':
            tim_read = [origin + datetime.timedelta(microseconds=dt) for dt in self.sp_data.coord("time").points]
        elif t_info.split(" ")[0] == 'milliseconds':
            tim_read = [origin + datetime.timedelta(milliseconds=dt) for dt in self.sp_data.coord("time").points]
        elif t_info.split(" ")[0] == 'minutes':
            tim_read = [origin + datetime.timedelta(minutes=dt) for dt in self.sp_data.coord("time").points]
        elif t_info.split(" ")[0] == 'hours':
            tim_read = [origin + datetime.timedelta(hours=dt) for dt in self.sp_data.coord("time").points]
        elif t_info.split(" ")[0] == 'weeks':
            tim_read = [origin + datetime.timedelta(weeks=dt) for dt in self.sp_data.coord("time").points]
        else:
            assert False, 'Wrong increment in coord("time")!'
        
        self.__tim_read__ = tim_read
            
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
        
        this_function = "overview"
        
        # TODO all the other plots
        
        list_of_plots=[]

        maxnumtemp = len(self.sp_data.coord("time").points)

        sp_masked_vals = self.sp_data.collapsed("time",iris.analysis.COUNT,function=lambda values: values.mask)
        sp_masked_vals.data.mask = sp_masked_vals.data==maxnumtemp
        sp_masked_vals.data = sp_masked_vals.data*0.+1.

        for d in self.dimensions:
            
            long_left_over = self.dimensions[self.dimensions!=d]
            short_left_over = np.array([sl[0:3] for sl in long_left_over])
                        
                        
            num_available_vals = self.sp_data.collapsed(d,iris.analysis.COUNT,function=lambda values: values>=np.min(self.sp_data.data.flatten()))
            if d in ["time"]:
                frac_available_vals = iris.analysis.maths.divide(num_available_vals,maxnumtemp)
            else:
                sp_agg_array=sp_masked_vals.collapsed(d,iris.analysis.SUM)   
                frac_available_vals = iris.analysis.maths.divide(num_available_vals,sp_agg_array.data)
            
            try:
                # plotting routine
                filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_frac_avail_" + "_".join(short_left_over) + "." + self.__output_type__
                list_of_plots.append(filename)
                x=Plot2D_2(frac_available_vals)
                
                fig = plt.figure()
                if "longitude" == d:     
                    gs = gridspec.GridSpec(1, 5)
                    ax = np.array([plt.subplot(gs[0, :-1]),plt.subplot(gs[0, -1])])
                    fig.set_figwidth(1.7*fig.get_figwidth())
                    fig.set_figheight(1.2*fig.get_figheight())
                elif "time" == d:
                    ax = [plt.subplot(1,1,1)]
                    fig.set_figheight(1.2*fig.get_figheight())
                elif "latitude" == d:
                    gs = gridspec.GridSpec(5, 1)
                    ax = np.array([plt.subplot(gs[:-1,0]),plt.subplot(gs[-1,0])])
                    fig.set_figheight(1.7*fig.get_figheight())
                x.plot(ax=ax, vminmax=[0.,1.], title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
                fig.savefig(filename)
                plt.close(fig.number)
                
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + ['DM_global', 'C3S_overview'],
                         str('Overview on ' + "/".join(long_left_over) + ' availablility of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                         '#C3S' + 'frav' + "".join(short_left_over) + self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)
            except:
                print('no figure done # 001')
            
        del sp_masked_vals
        del num_available_vals
        del frac_available_vals
        
        # histogram plot of available measurements
        all_data = self.sp_data.copy()
         
        try:
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
                     str('Full spatio-temporal histogram of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'histall' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
        except:
            print('no figure done #002')
        
        del all_data

        # dimension information
        lon_range = self.sp_data.coord("longitude").points
        lat_range = self.sp_data.coord("latitude").points
        tim_range = self.sp_data.coord("time").points
        t_info = str(self.sp_data.coord("time").units)
        
        lon_range_spec = utils.__minmeanmax__(lon_range)
        lat_range_spec = utils.__minmeanmax__(lat_range)
        tim_range_spec = utils.__minmeanmax__(tim_range)
        
        origin = datetime.datetime.strptime("_".join(t_info.split(" ")[-2:]), '%Y-%m-%d_%H:%M:%S')
        if t_info.split(" ")[0] == 'days':
            tim_range_spec_read = [origin + datetime.timedelta(days=dt) for dt in tim_range_spec]
        elif t_info.split(" ")[0] == 'seconds':
            tim_range_spec_read = [origin + datetime.timedelta(seconds=dt) for dt in tim_range_spec]
        elif t_info.split(" ")[0] == 'microseconds':
            tim_range_spec_read = [origin + datetime.timedelta(microseconds=dt) for dt in tim_range_spec]
        elif t_info.split(" ")[0] == 'milliseconds':
            tim_range_spec_read = [origin + datetime.timedelta(milliseconds=dt) for dt in tim_range_spec]
        elif t_info.split(" ")[0] == 'minutes':
            tim_range_spec_read = [origin + datetime.timedelta(minutes=dt) for dt in tim_range_spec]
        elif t_info.split(" ")[0] == 'hours':
            tim_range_spec_read = [origin + datetime.timedelta(hours=dt) for dt in tim_range_spec]
        elif t_info.split(" ")[0] == 'weeks':
            tim_range_spec_read = [origin + datetime.timedelta(weeks=dt) for dt in tim_range_spec]
        else:
            assert False, 'Wrong increment in coord("time")!'

        lon_freq = np.diff(lon_range)
        lat_freq = np.diff(lat_range)
        tim_freq = np.diff(tim_range)
        
        lon_freq_spec = utils.__minmeanmax__(lon_freq)
        lat_freq_spec = utils.__minmeanmax__(lat_freq)
        tim_freq_spec = utils.__minmeanmax__(tim_freq)
        
        overview_dict=collections.OrderedDict()
        
        overview_dict.update({'longitude range [' + str(self.sp_data.coord("longitude").units) + ']': collections.OrderedDict([("min",str(lon_range_spec[0])), ("max",str(lon_range_spec[2]))])})
        overview_dict.update({'longitude frequency [' + str(self.sp_data.coord("longitude").units) + ']': collections.OrderedDict([("min",str(lon_freq_spec[0])), ("average",str(lon_freq_spec[1])), ("max",str(lon_freq_spec[2]))])})
        overview_dict.update({'latitude range [' + str(self.sp_data.coord("latitude").units) + ']': collections.OrderedDict([("min",str(lat_range_spec[0])), ("max",str(lat_range_spec[2]))])})
        overview_dict.update({'latitude frequency [' + str(self.sp_data.coord("latitude").units) + ']': collections.OrderedDict([("min",str(lat_freq_spec[0])), ("average",str(lat_freq_spec[1])), ("max",str(lat_freq_spec[2]))])})
        overview_dict.update({'temporal range': collections.OrderedDict([("min",str(tim_range_spec_read[0])), ("max",str(tim_range_spec_read[2]))])})
        overview_dict.update({'temporal frequency [' + t_info.split(" ")[0] + ']': collections.OrderedDict([("min",str(tim_freq_spec[0])), ("average",str(round(tim_freq_spec[1],2))), ("max",str(tim_freq_spec[2]))])})
        
        # produce report
        self.__do_report__(content={"text":overview_dict,"plots":list_of_plots}, filename=this_function.upper())
        
        # save GCOS requirements
        self.__gcos_dict__.update({"Frequency":{"value":round(tim_freq_spec[1],2), "unit":t_info.split(" ")[0]}})
        self.__gcos_dict__.update({"Resolution":{"value":round(np.mean([lon_freq_spec[1],lat_freq_spec[1]]),2), "unit":str(self.sp_data.coord("longitude").units)}})
        
        # save values for other diagnostics
        self.__avg_timestep__ = tim_freq_spec[1]
        
        return
    
    def __do_report__(self,**kwargs):
        
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
    
    
    def __do_mean_var__(self):
        
        this_function = "mean & variability"
        
        list_of_plots = []
        
        maths = ["MEAN","STD_DEV","LOG_COV", "PERCENTILE", "CLIMATOLOGY"]
        
        percentiles = [1.,5.,10.,25.,1./3.*100,50.,2./3.*100,75.,90.,95.,99.]
        
        for d in self.dimensions:
            
            long_left_over = self.dimensions[self.dimensions!=d]
            short_left_over = np.array([sl[0:3] for sl in long_left_over])
        
            mean_std_cov=collections.OrderedDict()
            disp_min_max=collections.OrderedDict()
            disp_min_max.update({"abs_vals":np.array([np.nan])})
            disp_min_max.update({"diff_vals":np.array([np.nan])})

            for m in maths:
                
                if m == "PERCENTILE":
                    
                    try:
                        perc = self.sp_data.collapsed(d, iris.analysis.__dict__[m], percent=percentiles)
                        
                        for p in percentiles:
                            
                            loc_data = perc.extract(iris.Constraint(percentile_over_time=p))
                                                        
                            disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.nanpercentile(loc_data.data.data,[5,95])]),[0,100])})
                            
                            mean_std_cov.update({m + " " + str(int(round(p,0))) + " percent":loc_data})
                            
                            
                        del perc
                    except:
                        pass
                    
                elif m == "CLIMATOLOGY":
                    
                    try:
                        clim = self.sp_data.copy()
                        iris.coord_categorisation.add_month_number(clim, d, name='month_num')
                        clim_comp = clim.aggregated_by('month_num',iris.analysis.MEAN)
                        
                        clim_anom = clim.copy()
                        
                        del clim
                        
                        for mn in range(len(clim_anom.coord('month_num').points)):
                        
                            idx = clim_comp.coord('month_num').points.tolist().index(clim_anom.coord('month_num').points[mn])
                            clim_anom.data[mn,:,:] = clim_anom.data[mn,:,:]-clim_comp.data[idx,:,:]
                        
                        iris.coord_categorisation.add_year(clim_anom, d, name='year')
                        clim_anom = clim_anom.aggregated_by('year', iris.analysis.MEAN)
                        
                        for mon in clim_comp.coord('month_num').points:
                    
                            loc_data=clim_comp.extract(iris.Constraint(month_num=mon))
                            
                            mean_std_cov.update({m + " " + str(int(mon)):loc_data})
                        
                            disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.nanpercentile(loc_data.data.data,[5,95])]),[0,100])})
                        
                        del clim_comp
                        
                        for y in clim_anom.coord('year').points:
                            
                            loc_data=clim_anom.extract(iris.Constraint(year=y))
                            
                            mean_std_cov.update({"mean anomalies from " + m + " per year " + str(int(y)):loc_data})
                            
                            disp_min_max.update({"diff_vals":np.nanpercentile(np.concatenate([disp_min_max["diff_vals"],np.nanpercentile(loc_data.data.data,[5,95])]),[0,100])})
                                                        
                        del clim_anom
                            
                    except:
                        pass
                
                elif m == "MEAN":
                    
                    loc_data=self.sp_data.collapsed(d, iris.analysis.__dict__[m])
        
                    mean_std_cov.update({m:loc_data})
                    
                    disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.nanpercentile(loc_data.data.data,[5,95])]),[0,100])})
                    
                elif m == "STD_DEV":
        
                    mean_std_cov.update({m:self.sp_data.collapsed(d, iris.analysis.__dict__[m])})
                                 
                elif m == "LOG_COV": 
                    
                    mean_std_cov.update({m:iris.analysis.maths.log(iris.analysis.maths.divide(mean_std_cov["STD_DEV"],mean_std_cov["MEAN"]))})
         
                else:
                    raise ConfigurationError(m,"This maths functionality has to be implemented in _do_mean_var_.") 
            
            for m in mean_std_cov.keys():
                
                # plotting routine
                vminmax=None
                
                if np.any([m_typ in m for m_typ in ["MEAN", "PERCENTILE", "CLIMATOLOGY"]]):
                    vminmax=disp_min_max["abs_vals"]
                
                if np.any([m_typ in m for m_typ in ["anomalies"]]):
                    vminmax=disp_min_max["diff_vals"]
                    
                print d,m,vminmax
                
                if mean_std_cov[m] is not None:
                # this needs to be done due to an error in cartopy
                    filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "_".join(m.split(" ") )+ "_" + "_".join(short_left_over) + "." + self.__output_type__
                    list_of_plots.append(filename)
                    try:
                        x=Plot2D_2(mean_std_cov[m])
                        
                        fig = plt.figure()
                        if "longitude" == d:     
                            gs = gridspec.GridSpec(1, 5)
                            ax = np.array([plt.subplot(gs[0, :-1]),plt.subplot(gs[0, -1])])
                            fig.set_figwidth(1.7*fig.get_figwidth())
                            fig.set_figheight(1.2*fig.get_figheight())
                        elif "time" == d:
                            ax = [plt.subplot(1,1,1)]
                            fig.set_figheight(1.2*fig.get_figheight())
                        elif "latitude" == d:
                            gs = gridspec.GridSpec(5, 1)
                            ax = np.array([plt.subplot(gs[:-1,0]),plt.subplot(gs[-1,0])])
                            fig.set_figheight(1.7*fig.get_figheight())
                        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")",vminmax=vminmax)
                        fig.savefig(filename)
                        plt.close(fig.number)
                        # old plotting
    #                    fig = x.plot(summary_plot=True, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
    #                    fig.savefig(filename)
    #                    plt.close(fig)
                    
                        del mean_std_cov[m]
                        
                        ESMValMD("meta",
                                 filename,
                                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                                 str("/".join(long_left_over).title() + ' ' + m.lower() + ' values of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                                 '#C3S' + m + "".join(short_left_over) + self.__varname__,
                                 self.__infile__,
                                 self.diagname,
                                 self.authors)
                        
                    except Exception as e:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        print(exc_type, fname, exc_tb.tb_lineno)
                        print(mean_std_cov[m])
                        print('no figure done #003, blank figure')
                        
                        x=Plot2D_blank(mean_std_cov[m])
                        
                        fig = plt.figure()
                        if "longitude" == d:     
                            gs = gridspec.GridSpec(1, 5)
                            ax = np.array([plt.subplot(gs[0, :-1]),plt.subplot(gs[0, -1])])
                            fig.set_figwidth(1.7*fig.get_figwidth())
                            fig.set_figheight(1.2*fig.get_figheight())
                        elif "time" == d:
                            ax = [plt.subplot(1,1,1)]
                            fig.set_figheight(1.2*fig.get_figheight())
                        elif "latitude" == d:
                            gs = gridspec.GridSpec(5, 1)
                            ax = np.array([plt.subplot(gs[:-1,0]),plt.subplot(gs[-1,0])])
                            fig.set_figheight(1.7*fig.get_figheight())
                        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
                        fig.savefig(filename)
                        plt.close(fig.number)
                        # old plotting
    #                    fig = x.plot(summary_plot=True, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
    #                    fig.savefig(filename)
    #                    plt.close(fig)
                    
                        del mean_std_cov[m]
                        
                        ESMValMD("meta",
                                 filename,
                                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                                 str("/".join(long_left_over).title() + ' ' + m.lower() + ' values of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + '); Data can not be displayed due to cartopy error!'),
                                 '#C3S' + m + "".join(short_left_over) + self.__varname__,
                                 self.__infile__,
                                 self.diagname,
                                 self.authors)
                    
            del mean_std_cov
        
        # produce report
        self.__do_report__(content={"plots":list_of_plots}, filename=this_function.upper())

        return
    
    
    def __do_trends__(self):
        
        this_function = "trends & stability"
        
        list_of_plots = []
        
        # simple linear trend (slope) and p-values
        _,S,_,P = utils.__temporal_trend__(self.sp_data, pthres=1.01)        
        
        try:
            # plotting routines
            x=Plot2D_2(S)
            
            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_trend." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            ax = [plt.subplot(1,1,1)]
            fig.set_figheight(1.2*fig.get_figheight())
            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' slope values of ' + self.__varname__ + ' temporal trends per decade for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
            x=Plot2D_2(P)
    
            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_pvals." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            ax = [plt.subplot(1,1,1)]
            fig.set_figheight(1.2*fig.get_figheight())
            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' p-values for slopes of ' + self.__varname__ + ' temporal trends per decade for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
        except:
            print('no figure done #004')
        
        del P
        
        # linear trend (slope),breakpoints, and actual data after homogenization
        TempStab = utils.__TS_of_cube__(self.sp_data,
                                        dates=self.__tim_read__,
                                        breakpoint_method="CUMSUMADJ",
                                        max_num_periods=3,
                                        periods_method="autocorr",
                                        temporal_resolution=self.__avg_timestep__,
                                        min_avail_pts=2)
        
        # plotting routines
        try:
            # plotting the slope after break point correction
            x=Plot2D_2(TempStab["slope"])
    
            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_min_slope." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            ax = [plt.subplot(1,1,1)]
            fig.set_figheight(1.2*fig.get_figheight())
            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' slope values of ' + self.__varname__ + ' temporal trends per decade after breakpoint detection for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
            # plotting number of breakpoints
            x=Plot2D_2(TempStab["number_breakpts"])
            
            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_num_bps." + self.__output_type__
            list_of_plots.append(filename)
    
            fig = plt.figure()
            ax = [plt.subplot(1,1,1)]
            fig.set_figheight(1.2*fig.get_figheight())
            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' number of breakpoints of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
        except:
            print('no figure done #005')
        
        # plotting the version (2=deseason, 1=detrend, 0=neither, -1=not enough data available, -2=something went wrong)
#        x=Plot2D(TempStab["version"])
#        figS = x.plot(summary_plot=True, title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#
#        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_ver_trend." + self.__output_type__
#        list_of_plots.append(filename)
#        figS.savefig(filename)
#        plt.close(figS)
#        
#        ESMValMD("meta",
#                 filename,
#                 self.__basetags__ + ['DM_global', 'C3S_trend'],
#                 str("Latitude/Longitude" + ' versions of trend estimation of ' + self.__varname__ + ' temporal trends per decade after breakpoint detection for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                 '#C3S' + 'temptrend' + self.__varname__,
#                 self.__infile__,
#                 self.diagname,
#                 self.authors)
        try:
            x=Plot2D_2(TempStab["slope"]-S)
    
            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_hom_mean_diff." + self.__output_type__
            list_of_plots.append(filename)
    
            fig = plt.figure()
            ax = [plt.subplot(1,1,1)]
            fig.set_figheight(1.2*fig.get_figheight())
            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' slope difference after breakpoint reduction for ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
            x=Plot2D_2(cubestats.pearsonr(TempStab["homogenized"], self.sp_data, corr_coords="time"))
            
            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_hom_corr." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            ax = [plt.subplot(1,1,1)]
            fig.set_figheight(1.2*fig.get_figheight())
            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' correlation after breakpoint reduction for ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
        except:
            print('no figure done #006')
        
        del TempStab
        
        # produce report
        self.__do_report__(content={"plots":list_of_plots}, filename=this_function.upper())
        
        # update gcos
        self.__gcos_dict__.update({"Accuracy":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Stability":{"value":None, "unit":None}})
        
        return
    
    
    def __do_maturity_matrix__(self):
        
        this_function = "System maturity matrix"
        
        expected_input = self.__work_dir__ + os.sep + "smm_input" + os.sep + self.__basic_filename__ + "_smm_expert.csv"
        if not os.path.isfile(expected_input):
            try:
                os.makedirs(os.path.dirname(expected_input))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
            shutil.copy2(os.path.dirname(os.path.realpath(__file__)) + "/lib/predef/empty_smm_expert.csv", expected_input)
            print("************************************** WARNING **************************************")
            print("Expected " + this_function + " input file " + expected_input + " not found!")
            print("Created dummy file instead. Please fill in appropriate values and rerun!")
            print("(This won't fail if you do not do so and produce a white matrix!!!!)")
            print("************************************** WARNING **************************************")
            captionerror = True
        else:
            print("Processing " + this_function + " input file: " + expected_input) 
            captionerror = False

        # plotting routines
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_function.split()) + "." + self.__output_type__
        fig = do_smm_table(expected_input, os.path.dirname(os.path.realpath(__file__)) + "/lib/predef/smm_definitions.csv")
        fig.savefig(filename)
        plt.close(fig)
        
        if captionerror:
            caption = "The expected input file: " + expected_input + " was not found and an empty dummy file created, therefore this plot is blank. Please edit the expected file!"
        else:
            caption = str(this_function + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')


        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_SMM'],
                 caption,
                 '#C3S' + 'SMM' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # produce report
        self.__do_report__(content={"plots":[filename]}, filename="".join(this_function.upper().split()))
        
        return
    
    
    def __do_gcos_requirements__(self):
        
        this_function = "GCOS requirements"
        
        # plotting routines
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_function.split()) + "." + self.__output_type__
        fig = do_gcos_table(self.__varname__, self.__gcos_dict__, os.path.dirname(os.path.realpath(__file__)) + "/example_csvs/example_gcos_reference.csv")
        fig.savefig(filename)
        plt.close(fig)
        
        caption = str(this_function + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')
        
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_GCOS'],
                 caption,
                 '#C3S' + 'GCOS' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)
        
        # produce report
        self.__do_report__(content={"plots":[filename]}, filename="".join(this_function.upper().split()))
        
        return
