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

class AlbedoDiagnostic(BasicDiagnostics):
    """
    class to implement albedo diagnostics, like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """
    def __init__(self, **kwargs):
        super(AlbedoDiagnostic, self).__init__(**kwargs)
        
        self._project_info={}
        self._modtype = None 
        self._reftype = None
        self._plot_dir='.' + os.sep
        self._work_dir='.' + os.sep

        self._vartype = 'albedo' #default value as there must be one
        self.output_type = 'png'  #default ouput file type 
        self._changed=False

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(AlbedoDiagnostic, self).run_diagnostic()


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
        super(AlbedoDiagnostic, self).write_data()
        
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
                "AUX_Files_alb_QA4ECV"
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
                "AUX_Files_alb_QA4ECV"
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
        if proj_var == "alb":
            self._ref_data = self._load_cci_generic(self._ref_file, proj_var)
            self._ref_data.unit="-"
        else:
            assert False, 'Not supported yet'
            
            
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
                    with HiddenPrints():
                        TSA = TS(timearray, array, breakpoint_method="olssum",
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
        
        
        
######################### unused
    def _preprocess_observations(self, infile, mod, var,check_f = None,force=False):
        """
        preprocess observations to adapt to temporal and spatial resolution needed
        Parameters:
        -----------
        
        infile : string with path to file
        mod : model data; the prepocessing should mirror its specifications
        var : name of the variable within the file
        check_f : alternativeley check this folder and write infile.built.nc
        """
        
        #choose timestep and resolution 
        #TODO get this from mod_data or proj_info
        if self._field_type in ["T2Ms", "T3M"]:
            timestep="monthly"
        else:
            assert False, "unnkown field_type"
            
        resolution="T63"
        
        ofile=infile+'.built.nc'
        
        if (os.path.isfile(infile) or os.path.isfile(ofile)) and not force:
            data = self._load_cci_generic(ofile if os.path.isfile(ofile) else infile,var)
            return data
        elif (os.path.isfile(infile) or os.path.isfile(ofile)) and force:
            
            #adjust timestep
            thisfile = self._aggregate_timestep(ofile if os.path.isfile(ofile) else infile,timestep)
            
            #adjust resolution
            thisfile = self._aggregate_resolution(thisfile,resolution)
            
            subprocess.call(['cp',thisfile,ofile])
            
            os.remove(thisfile)
            
            data = self._load_cci_generic(ofile,var)
            return data
            
        elif not check_f is None:
            file_list, list_length = self._get_files_in_directory(check_f,'*/*/*/*.nc',False)
            
            #choose necessary files
            def loc_timestamp_split(str):
                return str.split('/')[-1].split('-')[0]
                
            file_timestamps = np.asarray(map(int,map(loc_timestamp_split,file_list)))
            ys=file_timestamps/10000000000
            low_date = min(ys)
            high_date = max(ys)
            
            curyear=low_date
                       
            file_list_agg=[]
            
            while not curyear == (high_date + 1):
                
                use = np.array(np.where(np.all([file_timestamps/10000000000 == curyear],0)))
                use=use[0][(np.argsort(file_timestamps[use[0]]))]
                file_list_cur=np.array(file_list)[use].tolist()
                
                if len(file_list_cur)>0:
                    #concatenate files
                    thisfile = self._aggregate_obs_from_files(file_list_cur)
                    
                    #adjust mask
                    thisfile = self._apply_sst_flags(thisfile)
                
                    #adjust timestep
                    thisfile = self._aggregate_timestep(thisfile,timestep)
                    
                    #adjust resolution
                    thisfile = self._aggregate_resolution(thisfile,resolution)
                    
                    #add to new file list
                    file_list_agg.append(thisfile)
                
                curyear = curyear + 1
                
                print(file_list_agg)

            thisfile = self._aggregate_obs_from_files(file_list_agg)
            
            [os.remove(tf) for tf in file_list_agg]

            subprocess.call(['cp',thisfile,ofile])
            
            os.remove(thisfile)
            
            data = self._load_cci_generic(ofile,var)
            return data
            
        else:
            print infile
            assert False, "cannot find any files!" 
            
    def _apply_sst_flags(self,infile,remove=True):
        """ apply flag to file """
        """ currenty only monthly """
        cdo=Cdo()
        oname=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmpname1=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmpname2=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmpname4=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1] #only if only sm is wanted
        cdo.selname("mask",input=infile,output=tmpname1,options='-L -f nc4 -b F32')
        cdo.selname("analysed_sst",input=infile,output=tmpname4,options='-L -f nc4 -b F32') #only if only sm is wanted
        cdo.setvrange(1,1,input=tmpname1,output=tmpname2,options='-L -f nc4 -b F32')
        os.remove(tmpname1)
        cdo.div(input=tmpname4 + " " + tmpname2,output=oname,options='-L -f nc4 -b F32') #only if only sm is wanted
        #cdo.div(input=infile + " " + tmpname3,output=oname,options='-L -f nc4 -b F32')
        os.remove(tmpname2)
        os.remove(tmpname4) #only if only sm is wanted
        
        if remove:
            os.remove(infile)
                
            
        return oname
        
###TODO CURRENTLY UNUSED

    def download_observations(self, year, force=True):
        """
        download observations and update local mirror with recent data

        The data is currently downloaded from the ICDC (http://icdc.zmaw.de/1.html)
        This requires a special account. An alternative approach would
        be to download the data directly from the SM CCI project
        and/or the CCI data protal

        TODO implement download from CCI or data portal
        """
        if not os.path.exists(self.raw_reference_directory):
            os.makedirs(self.raw_reference_directory)

        odir = self.raw_reference_directory + os.sep + year + os.sep
        if os.path.exists(odir):
            if force == False:
                print('No data download as data already available: ' + year)
                return

        if not os.path.exists(odir):
            os.makedirs(odir)
        else:  # empty directory
            os.system('rm -rf ' + odir + '*.nc')

        # download data via scp
        # if on the server, create only links, otherwise download using scp
        if os.path.exists(self.remote_dir):
            print('Linking data only as on server!')
            cmd = "ln -s " + self.remote_dir + year + '/*.nc' + ' ' + odir
            os.system(cmd)
        else:
            assert False, "not working yet with subprocess"
            cmd = self.user + '@' + self.server + ':' + self.remote_dir + year + '/*.nc' + ' ' + odir
            subprocess.call(["scp", "-r", cmd], shell=True)

    def update_observation(self, year, force=False):
        """
        download and update observations for a particular year

        Parameters
        ----------
        year : str
            year to process
        force : bool
            specifies if processing should be done in any case (true) or if
            the data is available, no processing will be made (false)
        """
        ofile_root = self.reference_preproc_directory + os.sep + 'CCI_SM_' + year

        self.download_observations(year, force=force)
        self.preprocess_observations(year, ofile_root, force=force)


