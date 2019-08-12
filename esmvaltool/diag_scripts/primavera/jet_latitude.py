import os
import logging

import matplotlib.pyplot as plt

import dask.array as da
import numpy as np

import iris
import iris.cube
import iris.analysis
import iris.util
import iris.coord_categorisation
import iris.quickplot as qp

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

logger = logging.getLogger(os.path.basename(__file__))


class JetLatitude(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.output_name = None
        self.target_grid = self.cfg.get('target_grid')
        self.grid_cube = None
        self.sftlf = None

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        lanczos_weight = np.load(
            '/home/users/panos/work_space/jetlat/python/LF_weights.npy'
        )
        for alias in data:
            ua = iris.load_cube(data[alias][0]['filename'])

            ua_filtered = np.apply_along_axis(
                lambda m: np.convolve(m, lanczos_weight, mode='same'),
                axis=ua.coord_dims('time')[0],
                arr=ua.core_data()
            )
            ua = ua.copy(ua_filtered)

            wind = ua.collapsed('latitude', iris.analysis.MAX)
            latitude = np.argmax(
                ua.data,
                axis=ua.coord_dims('latitude')[0]
            )
            del ua
            del ua_filtered
            latitude = wind.copy(latitude)
            latitude.var_name = 'lat'
            latitude.standard_name = 'latitude'
            latitude.long_name = 'Jet latitude'
            latitude.units = 'degrees_north'

            logger.debug(wind)
            logger.debug(latitude)

            self._compute_histogram(alias, wind)
            self._compute_histogram(alias, latitude)

    def _compute_histogram(self, alias, data):
        clim = data.aggregated_by('day_of_year', iris.analysis.MEAN)
        clim.remove_coord('time')
        clim.remove_coord('month_number')
        clim.remove_coord('day_of_month')
        clim.remove_coord('year')
        iris.util.promote_aux_coord_to_dim_coord(clim, 'day_of_year')
        logger.debug(clim)
        clim_fft = np.fft.rfft(clim.data)
        clim_fft[3:np.size(clim_fft)] = 0
        clim_fft = np.fft.irfft(clim_fft)
        clim = clim.copy(clim_fft)

        qp.plot(clim)
        plt.savefig(os.path.join(
            self.cfg[n.PLOT_DIR], '{}_{}_clim.png'.format(alias, data.var_name)))
        plt.close()

        anom = data - clim
        qp.plot(anom)
        plt.savefig(os.path.join(
            self.cfg[n.PLOT_DIR], '{}_{}_anom.png'.format(alias, data.var_name)))
        plt.close()


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        JetLatitude(config).compute()


if __name__ == '__main__':
    main()
