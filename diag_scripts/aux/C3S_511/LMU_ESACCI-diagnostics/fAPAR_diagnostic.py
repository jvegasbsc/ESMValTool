import os
import shutil
from netCDF4 import Dataset
from ESMValMD import ESMValMD
from diagnostic import BasicDiagnostics
from TempStab.TempStab import TempStab as TS
import matplotlib.pyplot as plt
import numpy as np
from geoval.core.mapping import SingleMap

import sys

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout
        
class fAPARDiagnostic(BasicDiagnostics):
    """
    class to implement albedo diagnostics, like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """
    def __init__(self, **kwargs):
        super(fAPARDiagnostic, self).__init__(**kwargs)
        
        self._project_info={}
        self._modtype = None 
        self._reftype = None
        self._plot_dir='.' + os.sep
        self._work_dir='.' + os.sep

        self._vartype = 'fAPAR' #default value as there must be one
        self.output_type = 'png'  #default ouput file type 
        self._changed=False

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(fAPARDiagnostic, self).run_diagnostic()

    def _specific_diag(self, percentile=True):
        """
        Diagnostic management
        """
        
        if True:
            self._temporal_stability()


    def write_data(self,plot=True):
        """
        write data
        """
        super(fAPARDiagnostic, self).write_data()
        
        if '_ref_min_trend' and '_mod_min_trend'and \
            '_ref_num_bp' and '_mod_num_bp' in self.__dict__.keys():
            self._plot_TSA_maps()
        else:
            print 'No temporal stability to plot!'
        
        
    def _load_model_data(self):
        """ load albedo model data """
        
        mod_info = Dataset(self._mod_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        if not (lat == self.reso["lat"] and lon == self.reso["lon"]):

            grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                str(self.reso["lon"])]

            newfile = self._mod_file + "." + grid + "built.nc"
            newfile = newfile.split("/")
            newdir = (self._work_dir if self._work_dir[-1] ==
                      os.sep else self._work_dir + os.sep) +\
                "AUX_Files_fAPAR_QA4ECV"
            newfile = newdir + os.sep + newfile[-1]

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._mod_file, grid, remove=False)
                if not os.path.exists(newdir):
                    os.makedirs(newdir)
                shutil.copy2(tempfile, newfile)
                os.remove(tempfile)

            self._mod_file = newfile

        # load data
        proj_var = self._project_info['RUNTIME']['currDiag'].get_variables()[0]
        self._mod_data = self._load_cmip_generic(self._mod_file, proj_var)
        self._mod_data.unit="-"
        
        
    def _load_observation_data(self):
        """ load obs data """
        
        mod_info = Dataset(self._ref_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        reso = {"lat": lat, "lon": lon}

        if str(lat) + "-" + str(lon) not in self.resoLUT.keys():

            grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                str(self.reso["lon"])]

            newfile = self._ref_file + "." + grid + "built.nc"
            newfile = newfile.split("/")
            newdir = (self._work_dir if self._work_dir[-1] ==
                      os.sep else self._work_dir + os.sep) +\
                "AUX_Files_fAPAR_QA4ECV"
            newfile = newdir + os.sep + newfile[-1]

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._ref_file, grid, remove=False)
                if not os.path.exists(newdir):
                    os.makedirs(newdir)
                shutil.copy2(tempfile, newfile)
                os.remove(tempfile)

            self._ref_file = newfile

            reso = self.reso

        self.reso = reso

        proj_var = self._project_info['RUNTIME']['currDiag'].get_variables()[0]
        if proj_var == "fAPAR":
            self._ref_data = self._load_cci_generic(self._ref_file, proj_var)
            self._ref_data.unit="-"
        else:
            assert False, 'Not supported yet'


    def _load_regionalization_shape(self):
        """ load shape data """
        self._reg_shape = self._load_shape_generic(self._reg_file)
        
        
    def _temporal_stability(self):
        """
        calculate temporal stability for model and observational dataset
        (and compare these)
        """

        print('   temporal stability analysis ...')
        
        def loc_TSA_fun(array,**kwargs):
            
            RES = None
            count=0
            
            timearray = kwargs.get('dates', None)
            
            if timearray is not None:
                try:
#                    with HiddenPrints():
                    TSA = TS(timearray, array, breakpoint_method="CUSUMADJ",
                             deseason=True, max_num_periods=3,
                             periods_method="autocorr",
                             temporal_resolution = "monthly")
                    RES = TSA.analysis(homogenize=True)
                except:
                    count += 1
            else:
                "Error in timearray."
                
            if RES is not None:
                slope_diff = (np.log(abs(RES["homogenized_trend"]["slope"]))-
                              np.log(abs(RES["original_trend"]["slope"]))) / \
                             np.log(abs(RES["original_trend"]["slope"])) * \
                             np.sign(RES["homogenized_trend"]["slope"]) * \
                             np.sign(RES["original_trend"]["slope"])
                return np.atleast_1d(np.array([slope_diff*100.,
                                               len(RES["breakpoints"]),
                                               count]))
            
            else: 
                return np.atleast_1d(np.array([np.nan, np.nan, count]))
        
        self._ref_min_trend = self._ref_data.get_percentile(0.)
        self._ref_num_bp = self._ref_min_trend.copy()
        self._mod_min_trend = self._mod_data.get_percentile(0.)
        self._mod_num_bp = self._mod_min_trend.copy()
        
        res_ref = np.apply_along_axis(loc_TSA_fun, 0, self._ref_data.data,
                                      dates = self._ref_data.date)
        print("    Temproal stability analysis not applicable for " + 
              str(int(np.sum(res_ref[2,:,:]) -
                  np.sum(self._ref_min_trend.data.mask))) + " of " +
              str(int(np.sum(np.logical_not(self._ref_min_trend.data.mask)))) +
              " unmasked time series in reference data!")
        
        res_mod = np.apply_along_axis(loc_TSA_fun, 0, self._mod_data.data,
                                      dates = self._mod_data.date)
        print("    Temproal stability analysis not applicable for " + 
              str(int(np.sum(res_mod[2,:,:]) - 
                  np.sum(self._mod_min_trend.data.mask))) + " of " +
              str(int(np.sum(np.logical_not(self._mod_min_trend.data.mask)))) +
              " unmasked time series in model data!")
        
        ref_mask = np.isnan(res_ref[0,:,:]) + self._ref_min_trend.data.mask
        self._ref_min_trend.data = np.ma.array(data=res_ref[0,:,:],
                                               mask=ref_mask)
        self._ref_min_trend.unit = "%"
        
        mod_mask = np.isnan(res_mod[0,:,:]) + self._mod_min_trend.data.mask
        self._mod_min_trend.data = np.ma.array(data=res_mod[0,:,:],
                                               mask=mod_mask)
        self._mod_min_trend.unit = "%"
        
        self._ref_num_bp.data = np.ma.array(data=res_ref[1,:,:],
                                            mask=ref_mask)
        self._mod_num_bp.data = np.ma.array(data=res_mod[1,:,:],
                                            mask=mod_mask)
        
    def _plot_TSA_maps(self):
        """
        plot temporal stability information
        """
        f = plt.figure(figsize=(20, 14))
        ax1 = f.add_subplot(221)
        ax2 = f.add_subplot(222)
        ax3 = f.add_subplot(223)
        ax4 = f.add_subplot(224)

        def submap(data, ax, title, vmin, vmax, cmap,
                   ctick={'ticks': None, 'labels': None}):
            Map = SingleMap(data,
                            backend=self.plot_backend,
                            show_statistic=True,
                            savefile=None,
                            ax=ax,
                            show_unit=True)
            Map.plot(title=title,
                     show_zonal=False,
                     show_histogram=False,
                     show_timeseries=False,
                     nclasses=self.plot_nclasses,
                     colorbar_orientation=self.plot_cborientation,
                     show_colorbar=self.plot_cbshow,
                     cmap=cmap,
                     vmin=vmin,
                     vmax=vmax,
                     proj_prop=self.cfg.projection,
                     ctick_prop=ctick,
                     drawparallels=True,
                     titlefontsize=self.plot_tfont)
            
        # TODO should be slope diff
        submap(self._ref_min_trend, ax=ax1, title="relative log slope change to minimal detectable trend (" + self.refname + ")", vmin=-25, vmax=+25, cmap='RdBu')
        submap(self._mod_min_trend, ax=ax2, title="relative log slope change to minimal detectable trend (" + self.modname + ")", vmin=-25, vmax=+25, cmap='RdBu')
        
        lab_max = np.max([np.nanmax(self._ref_num_bp.data), np.nanmax(self._mod_num_bp.data)])
        tens = 10**np.floor(np.log10(lab_max))
        lab_max = np.ceil(lab_max/tens)*tens
        
        submap(self._ref_num_bp, ax=ax3, title="number of detected breakpoints (" + self.refname + ")", vmin=0, vmax=lab_max, cmap='summer')
        submap(self._mod_num_bp, ax=ax4, title="number of detected breakpoints (" + self.modname + ")", vmin=0, vmax=lab_max, cmap='summer')

        titlesup = "temporal stability of " + \
            self.refname + " and " + \
            self.modname
        oname = self._get_output_rootname() + \
            "_temporal_stability" + '.' + self.output_type

        f.suptitle(titlesup)

        if os.path.exists(oname):
            os.remove(oname)
        f.savefig(oname, dpi=self.plot_dpi)

        plt.close(f.number)  # close figure for memory reasons!
        del f

# TODO!
        ESMValMD("both",
                 oname,
                 self._basetags + ['DM_global', 'PT_geo', 'ST_corr',
                                   self.refname, self.modname],
                 str('Pixelwise correlation of temporal trend between ' +
                     self.modname + ' and ' + self.refname + ' ' +
                     self._vartype + (' anomalies' if 'anomalytrend' in
                                      self.cfg.__dict__.keys() and
                                      self.cfg.anomalytrend else '') +
                     ' (left). The right panel describes ' +
                     'the pattern of corresponding p-values.'),
                 '#ID' + 'TSA' + self.var,
                 ",".join(self._infiles))
        