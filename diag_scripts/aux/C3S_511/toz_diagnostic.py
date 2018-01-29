import os
import shutil
from netCDF4 import Dataset
from diagnostic import BasicDiagnostics


class TOZDiagnostic(BasicDiagnostics):
    """
    class to implement TOZ diagnostics,
    like e.g. global means, global differences, RMSD etc.

    """

    def __init__(self, **kwargs):
        super(TOZDiagnostic, self).__init__(**kwargs)

        self._project_info = {}
        self._modtype = None
        self._reftype = None
        self._plot_dir = '.' + os.sep
        self._work_dir = '.' + os.sep

        self._vartype = 'TOZ'  # default value as there must be one
        self.output_type = 'png'  # default ouput file type
        self._changed = False

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(TOZDiagnostic, self).run_diagnostic()

    def _specific_diag(self):
        """
        Diagnostic management
        """

    def write_data(self, plot=True):
        """
        write data
        """
        super(TOZDiagnostic, self).write_data()

	def _load_model_data(self):
		""" load model data """

		print("test in dummy_diagnostics.py, load_model_data")
		
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
                "AUX_Files_TOZ"
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
                "AUX_Files_TOZ"
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
        if proj_var == "toz":
            self._ref_data = self._load_cci_generic(self._ref_file, proj_var)
        else:
            assert False, 'Not supported yet'


	def load_data(self):
		"""
		load model and observation data
		actual implementation needs to be part of child class
		"""
		
		print("finally in load_data")
		
		self._write_loading_header()

		if not self._changed:
			self._start_time, self._stop_time = None, None

		self._load_observation_data()
		self._load_model_data()

		# compatible time range
		# if self._project_info['RUNTIME']['currDiag'].get_variables()[0] in
		# ['natural_grass', 'managed_grass_and_crops', 'bare_soil', 'shrubs',
		# 'forest','lc'] else "month")
		self._start_time, self._stop_time, self._changed = \
			self._adjust_time_range(scale="month")

		if 'start_year' in self.cfg.__dict__.keys():
			if self.cfg.start_year > self._start_time.year:
				self._start_time = self._start_time.replace(
					year=self.cfg.start_year, month=01, day=01)
				self._changed = True
		if 'stop_year' in self.cfg.__dict__.keys():
			if self.cfg.stop_year < self._stop_time.year:
				self._stop_time = self._stop_time.replace(
					year=self.cfg.stop_year, month=12, day=31)
				self._changed = True

		if self._changed and not (self.var in ["baresoilFrac",
                                               "grassNcropFrac",
                                               "shrubNtreeFrac"]):
			# reorganize data
			self._mod_data.apply_temporal_subsetting(
				self._start_time, self._stop_time)
			self._ref_data.apply_temporal_subsetting(
				self._start_time, self._stop_time)

		# compatible spatial masks
		self._check_mask()
		self._ts = [x.date() for x in self._mod_data.date]

		if "climatologies" in self.cfg.__dict__.keys() and \
				self.cfg.climatologies:
			self._mts = list(set([int(x.month) for x in self._ts]))
			if len(self._mts) is not 12:
				self.cfg.climatologies = False
				print('   time series to short; ' +
                      'no calculation of climatologies ...')

		if 'write_preprocessed_data' in self.cfg.__dict__.keys():
			if self.cfg.write_preprocessed_data:
				self._save_p_data()

		self._ts = [x.date() for x in self._mod_data.date]

#        if 'write_preprocessed_data' in self.cfg.__dict__.keys():
#            if self.cfg.write_preprocessed_data:
#                self._save_p_data()

		# make sure that only used netcdf files set an attribute
		# ending on "_file"
		self._infiles = [getattr(self, s) for s in self.__dict__.keys()
                         if ("_mod_file" in s or "_ref_file" in s)]
		self._allfiles = [getattr(self, s) for s in self.__dict__.keys()
                          if ("_file" in s)]

	def calc_mean(self):
		
		print("this is a test to see if I can figure thingsout")