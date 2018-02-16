#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:17:48 2018

"""

import imp
import numpy as np

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
        