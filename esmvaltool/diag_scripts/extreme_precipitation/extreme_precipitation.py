# =====================================
# generalised_extreme_value_analysis.py
# =====================================
#
# Apply GEV analysis to daily precipitation output by PRIMAVERA
# Stream 1 model simulations, grid point-by-grid point (WP1).
# Additional script available to aggregate GEV results over
# large river basins within Europe (WP10).
#
# Notes:
#   * Parametric 1-day block maxima method applied seasonally.
#   * GEV analysis run within R interface ('extRemes').
#   * EC-Earth model grid coords are transformed herein.
#   * CNRM-CERFACS model is transformed to regular grid and
#     transformed data are output herein.
#   * CMCC model data are preprocessed from 6-hourly output
#     using 'preprocess_cmcc_day_pr.py'.
#
# Alexander J. Baker, UREAD
# 19/10/2017
# alexander.baker@reading.ac.uk



# ----------
# USER INPUT
# ----------
r_period = [2.,5.,10.,20.,30.,50.]  # Return periods (in years)
sig = 0.05                          # Confidence interval



#############################################################


import os
import logging
import string

import matplotlib
matplotlib.use('Agg')  # noqa
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.font_manager
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox

from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont

import iris
import iris.cube
import iris.analysis
import iris.util

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
import numpy as np

# R interface
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri

extRemes = importr('extRemes')
rpy2.robjects.numpy2ri.activate()
R['options'](warn = -1)
R('rm(list = ls())')


logger = logging.getLogger(os.path.basename(__file__))


class ExtremePrecipitation(object):

    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.return_period = self.cfg['return_period']
        self.confidence_interval = self.cfg['confidence_interval']

        # --------------
        # GEV PARAMETERS
        # --------------
        self.gev_par_sym = ['mu', 'sigma', 'xi']
        self.gev_par_name = dict(mu='location', sigma='scale', xi='shape')
        self.r_period_name = np.array(self.return_period).astype(int).astype(str)
        for r in range(len(self.r_period_name)):
            self.r_period_name[r] = '{}-year level'.format(self.return_period[r])
        self.return_period = rpy2.robjects.IntVector(self.return_period)

    def compute(self):
        for filename in self.filenames:
            logger.info('Processing %s', filename)
            cube = iris.load_cube(filename)
            logger.info(cube)
            cube.data
            logger.info('GEV analysis...')

            fevd = dict()
            for par in self.gev_par_sym:
                fevd[par] = iris.cube.CubeList()

            rl = dict()
            for r in self.r_period_name:
                rl[r] = iris.cube.CubeList()

            lats = cube.coord('latitude').points
            lons = cube.coord('longitude').points
            logger.info(cube)

            for slice_cube in cube.slices_over(('latitude', 'longitude')):
                if np.any(slice_cube.data[:]):
                    bmax_ll = R.matrix(slice_cube.data)
                    evdf = extRemes.fevd(bmax_ll, units=slice_cube.units.origin)
                    if evdf.rx2('results').rx2('par').rx2('location')[0] > 0. and\
                        evdf.rx2('results').rx2('par').rx2('scale')[0] > 0.: # -ve mu/sigma invalid
                        for par in self.gev_par_sym:
                            val_ll = evdf.rx2('results').rx2('par').rx2(self.gev_par_name[par])[0]
                            fevd[par] = iris.cube.Cube(
                                val_ll,
                                long_name=par,
                                units=slice_cube.units,
                                aux_coords_and_dims=(
                                    (cube.coord('latitude'), None),
                                    (cube.coord('longitude'), None),
                                    (cube.coord('season'), None)
                                ),
                            )
                        r_level = extRemes.return_level(evdf, return_period=self.return_period,
                                                        qcov=extRemes.make_qcov(evdf))
                        for r in range(len(r_level)):
                            rl[self.r_period_name[r]] = iris.cube.Cube(
                                r_level[r],
                                long_name=self.r_period_name[r],
                                units=slice_cube.units,
                                aux_coords_and_dims=(
                                    (cube.coord('latitude'), None),
                                    (cube.coord('longitude'), None),
                                    (cube.coord('season'), None)
                                ),
                            )
                else:
                    for par in self.gev_par_sym:
                        fevd[par] = iris.cube.Cube(
                            np.nan,
                            long_name=par,
                            units=slice_cube.units,
                            aux_coords_and_dims=(
                                (cube.coord('latitude'), None),
                                (cube.coord('longitude'), None),
                                (cube.coord('season'), None)
                            ),
                        )
                    for r in self.self.r_period_name:
                        rl[r] = iris.cube.Cube(
                            np.nan,
                            long_name=r,
                            units=slice_cube.units,
                            aux_coords_and_dims=(
                                (cube.coord('latitude'), None),
                                (cube.coord('longitude'), None),
                                (cube.coord('season'), None)
                            ),
                        )

            # Output results
            results_subdir = os.path.join(results_dir, PROJECT,OUT,inst_n,
                                          ENSEMBLE[inst][model],forcing,runID)
            if not os.path.exists(results_subdir):
                os.makedirs(results_subdir)
            for par in self.gev_par_sym:
                par_cube = fevd[par].merge_cube()
                par_fn = '_'.join([ENSEMBLE[inst][model],forcing,runID,season,par,domain])
                par_ffp = os.path.join(results_subdir,par_fn+frmt)
                iris.save(par_cube,par_ffp)
            rl_cubelist = iris.cube.CubeList()
            for r in range(len(r_level)):
                rl_cubelist.append(iris.cube.Cube(rl[r_period_name[r]],
                    long_name = r_period_name[r],units = cube.units.origin,
                    dim_coords_and_dims = [(cube.coord('latitude'),0),
                                           (cube.coord('longitude'),1)]))
            rl_fn = '_'.join([ENSEMBLE[inst][model],forcing,runID,season,'rl',domain])
            rl_ffp = os.path.join(self.cfg[n.WORK_DIR], rl_fn+frmt)
            print(' saving '+rl_fn)
            iris.save(rl_cubelist,rl_ffp)



def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        ExtremePrecipitation(config).compute()


if __name__ == "__main__":
    main()