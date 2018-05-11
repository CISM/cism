#!/usr/bin/env python

# This script plots the grounding line position for all 3 MISMIP3d experiments at the end of each experiment.
# This script requires the user to have run the python script "mismip3dWriteGL.py".


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

####################################
# Function used later in the code #
####################################


def glplot(ncfile, times, colora, label):
    """
    add a plot of grounding line points to current axes.
    makes use of the numpy.ma.MaskedArray when reading xGL,yGL
    """
    ncid  = Dataset(ncfile, 'r')
    time  = ncid.variables["time"][:]
    lxmax = 0.0
    lxmin = 800.0
    for i in range(0, len(times)):
        seq   = (time == times[i])
        xGL   = ncid.variables["xGL"][:, seq]*1e-3
        lxmax = max(np.max(xGL), lxmax)
        lxmin = min(np.min(xGL), lxmin)
        yGL   = ncid.variables["yGL"][:, seq]*1e-3
        plt.plot(xGL, yGL, 's', ms=3, mfc=colora[i],
                 mec=colora[i], label=label + ', t = ' + format(times[i]))
    return lxmin, lxmax


########
# Code #
########

model = '_cism'


plt.figure(figsize=(7, 7))

fileStd    = 'Stnd/Stnd' + model + '.nc'
ncidStd    = Dataset(fileStd,'r')
timeStd    = ncidStd.variables["time"][:]
xmin, xmax = glplot(fileStd, [timeStd[-1]], ['black'], 'Stnd')
xminplot   = xmin

fileP75S   = 'P75S/P75S' + model + '.nc'
ncidP75S   = Dataset(fileP75S,'r')
timeP75S   = ncidP75S.variables["time"][:]
xmin, xmax = glplot(fileP75S, [timeP75S[-1]], ['red'], 'P75S')
xminplot   = xmin
xmaxplot   = xmax

fileP75R   = 'P75R/P75R' + model + '.nc'
ncidP75R   = Dataset(fileP75R,'r')
timeP75R   = ncidP75R.variables["time"][:]
plt.xlim([xminplot-50.0, xmaxplot+50.0])
xmin, xmax = glplot(fileP75R, [timeP75R[-1]], ['blue'], 'P75R')

plt.legend(frameon=True, borderaxespad=0, loc='right')
plt.xlabel(r'$x$ (km)')
plt.ylabel(r'$y$ (km)')

# Saving the figure.
plt.savefig("mismip3dPlotGL.pdf")
