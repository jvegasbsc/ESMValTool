#/usr/bin/env python2
# -*- coding: utf-8 -*-

import os

import iris
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import plot

# Settings
NUMBER_OF_CUBES = 4
DIR = os.path.expanduser('~/C3S_tests')
LATLON_CUBE_PATH = os.path.join('test_datasets', 'latlon.nc')
LATTIME_CUBE_PATH = os.path.join('test_datasets', 'lattime.nc')
LONTIME_CUBE_PATH = os.path.join('test_datasets', 'lontime.nc')
LATLON_OUT_PATH = os.path.join(DIR, 'latlon_plot_')
LATTIME_OUT_PATH = os.path.join(DIR, 'lattime_plot')
LONTIME_OUT_PATH = os.path.join(DIR, 'lontime_plot')
PLOT_TYPE = 'pdf'
if not os.path.exists(DIR):
    os.makedirs(DIR)

# Latlon cubes
LATLON_CUBES = []
LATLON_CUBE = iris.load_cube(LATLON_CUBE_PATH)
LCUBE = [LATLON_CUBE]
LATLON_CUBES.append(LATLON_CUBE)
for factor in range(1, NUMBER_OF_CUBES + 1):
    LATLON_CUBES.append(LCUBE * factor)

# Time cube
LATTIME_CUBE = iris.load_cube(LATTIME_CUBE_PATH)
LONTIME_CUBE = iris.load_cube(LONTIME_CUBE_PATH)

# Matplotlib
gs = gridspec.GridSpec(1, 5)

# Lattime tests
FIG = plt.figure()
AX = [plt.subplot(gs[0, :-1]), plt.subplot(gs[0, -1])]
plot2d = plot.Plot2D(LATTIME_CUBE)
plot2d.plot(ax=AX, summary_plot=True)
filepath = LATTIME_OUT_PATH + '.' + PLOT_TYPE
FIG.savefig(filepath)
print("*** Wrote " + filepath + " ***")
plt.close

# Latlon tests
for (idx, cubes) in enumerate(LATLON_CUBES):
    FIG = plt.figure()
    AX = plt.subplot()
    print(idx)
    plot2d = plot.Plot2D(cubes)
    plot2d.plot(ax=AX)
    filepath = LATLON_OUT_PATH + str(idx) + '.' + PLOT_TYPE
    FIG.savefig(filepath)
    print("*** Wrote " + filepath + " ***")
    plt.close()

# Lontime tests
FIG = plt.figure()
AX = [plt.subplot(gs[0, :-1]), plt.subplot(gs[0, -1])]
plot2d = plot.Plot2D(LONTIME_CUBE)
plot2d.plot(ax=AX, summary_plot=True)
filepath = LONTIME_OUT_PATH + '.' + PLOT_TYPE
FIG.savefig(filepath)
print("*** Wrote " + filepath + " ***")
plt.close()
