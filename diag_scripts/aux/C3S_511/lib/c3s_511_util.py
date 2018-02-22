#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:17:48 2018

"""

import imp
import numpy as np
from scipy import stats
import cf_units
import matplotlib.pyplot as plt


def __getInfoFromFile__(filename):
    """
    routine to read cfg file
    """
    f = open(filename)
    __config__ = imp.load_source('cfg', '', f)
    f.close()
    return __config__


def  __minmeanmax__(array):
    """
    calculate minimum, maximum and average of array
    """
    return (np.min(array),np.mean(array),np.max(array))


def __temporal_trend__(cube, pthres=1.01):
    """
    calculate temporal trend of the data over time
    the slope of the temporal trend has unit [dataunit/day]
    Parameters
    ----------
    return_object : bool
        specifies if a C{Data} object shall be returned [True]
        or if a numpy array shall be returned [False]
    pthres : float
        specifies significance threshold; all values above this threshold will be masked
    Returns
    -------
    The following variables are returned:
    correlation, slope, intercept, p-value
    """
    # time difference in days
    dt = np.asarray(cube.coord('time').points - cube.coord('time').points[0]
                    ).astype('float')/365.2425/10
    # ensure that a masked array is used
    x = np.ma.array(dt, mask=dt != dt)

    R, S, I, P, C = __corr_single__(cube, x, pthres=pthres)

    R.long_name = cube.long_name + ' (correlation)'
    S.long_name = cube.long_name + ' ($\partial x / \partial t$)'
    I.long_name = cube.long_name + ' (offset)'
    P.long_name = cube.long_name + ' (p-value)'

    if S.units in [None,'no_unit','1','unknown']:
        S.units = cf_units.Unit('0.1 year-1')
    else:
        S.units += cf_units.Unit(str(S.units) + ' 0.1 year-1')
    R.units = '-'
    I.units = cube.units
    P.units = '-'
    return R, S, I, P
    

def __corr_single__(cube, x, pthres=1.01):
    """
    The routine correlates a data vector with all data of the
    current object. It returns several correlation measures
    Parameters
    ----------
    cube : 
    x : ndarray
        the data vector correlations should be calculated with,
        numpy array [time]
    pthres : float
        significance threshold. All values below this
        threshold will be returned as valid
    """

    if cube.ndim != 3:
        raise ValueError('Invalid geometry!')

    nt, ny, nx = cube.shape

    if nt != len(x):
        raise ValueError('Inconsistent geometries')

    # check if 'x' is a masked array where the mask has the same
    # size as the data. If this is not the case, then set it
    if isinstance(x, np.ma.core.MaskedArray):
        if isinstance(x.mask, np.ndarray):
            if x.mask.shape != x.data.shape:
                raise ValueError('Invalid mask geometry!')
        else:
            raise ValueError(
                'The mask of the dataset needs to be an array!')
    else:
        raise ValueError('Expect masked array as input in corr_single')

    # get data with at least one valid value
    dat, msk = __get_valid_data__(cube, mode='thres', thres=3)
    xx, n = dat.shape
    print('   Number of grid points: ', n)

    R = np.ones((ny, nx)) * np.nan  # output matrix for correlation
    P = np.ones((ny, nx)) * np.nan  # output matrix for p-value
    S = np.ones((ny, nx)) * np.nan  # output matrix for slope
    I = np.ones((ny, nx)) * np.nan  # output matrix for intercept
    CO = np.ones((ny, nx)) * np.nan  # output matrix for covariance

    R.shape = (-1)
    S.shape = (-1)
    P.shape = (-1)
    I.shape = (-1)
    CO.shape = (-1)

    print 'Calculating correlation ...'
    res = [stats.mstats.linregress(x, dat[:, i]) for i in xrange(n)]

    res = np.asarray(res)
    slope = res[:, 0]
    intercept = res[:, 1]
    r_value = res[:, 2]
    p_value = res[:, 3]
#    std_err = res[:, 4]

    R[msk] = r_value
    P[msk] = p_value
    I[msk] = intercept
    S[msk] = slope
    R.shape = (ny, nx)
    P.shape = (ny, nx)
    I.shape = (ny, nx)
    S.shape = (ny, nx)

    #--- prepare output data objects
    Rout = cube[0,:,:].copy()  # copy object to get coordinates
    Rout.long_name = 'correlation'
    msk = (P > pthres) | (np.isnan(R))
    #msk = np.zeros_like(R).astype('bool')
    Rout.data = np.ma.array(R, mask=msk).copy()
    Rout.units = None

    Sout = cube[0,:,:].copy()  # copy object to get coordinates
    Sout.long_name = 'slope'
    Sout.data = np.ma.array(S, mask=msk).copy()
    Sout.units = cube.units

    Iout = cube[0,:,:].copy()  # copy object to get coordinates
    Iout.long_name = 'intercept'
    Iout.data = np.ma.array(I, mask=msk).copy()
    Iout.units = cube.units

    Pout = cube[0,:,:].copy()  # copy object to get coordinates
    Pout.long_name = 'p-value'
    Pout.data = np.ma.array(P, mask=msk).copy()
    Pout.units = None

    Cout = cube[0,:,:].copy()  # copy object to get coordinates
    Cout.long_name = 'covariance'
    # currently not supported: covariance!
    Cout.data = np.ma.array(np.ones(P.shape) * np.nan, mask=msk).copy()
    Cout.units = None

    return Rout, Sout, Iout, Pout, Cout
    

def __get_valid_data__(cube, mode='all', thres=-99):
    """
    this routine calculates from the masked array
    only the valid data and returns it together with its
    coordinate as vector
    valid means that ALL timestamps need to be valid!
    Parameters
    ----------
    return_mask : bool
        specifies if the mask applied to the original data should
        be returned as well
    mode : str
        analysis mode ['all','one','thres']
        'all': all timestamps need to be valid
        'one': at least a single dataset needs to be valid
        'thres' : number of valid timesteps needs to be abovt a threshold
    thres : int
        threshold for minimum number of valid values (needed when mode=='thres')
    """

    if mode == 'thres':
        assert thres > 0, 'Threshold needs to be > 0!'

    if cube.ndim == 3:

        n = len(cube.coord("time").points)

        # vectorize the data

        data = cube.data.reshape(n, -1)
        # set pixels with NaN to invalid
        data.mask[np.isnan(data.data)] = True

        # extract only valid (not masked data)
        if mode == 'all':
            # identify all ngrid cells where all timesteps are valid
            msk = np.sum(~data.mask, axis=0) == n
        elif mode == 'one':
            # identify ONE grid cell where all timesteps are valid
            msk = np.sum(~data.mask, axis=0) > 0
        elif mode == 'thres':
            msk = np.sum(~data.mask, axis=0) > thres
        else:
            raise ValueError('Invalid option in get_valid_data() %s' %
                             mode)

        data = data[:, msk]

    elif cube.ndim == 2:
        data = cube.data.reshape(-1)
        msk = ~data.mask
        data = data[msk]
    else:
        raise ValueError('Unsupported dimension!')

    return data, msk



    