#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:17:48 2018

"""

import imp

def __getInfoFromFile__(filename):
    """
    routine to read cfg file
    """
    f = open(filename)
    __config__ = imp.load_source('cfg', '', f)
    f.close()
    return __config__
        