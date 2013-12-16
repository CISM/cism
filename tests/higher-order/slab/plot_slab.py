#!/usr/bin/env python
# This script plots result of an experiment with an ice "slab" on an inclined plane.
# Written by Matt Hoffman, Dec. 16, 2013
# Test case described in sections 5.1-2 of:
# J.K. Dukoqicz, 2012. Reformulating the full-Stokes ice sheet model for a more efficient computational solution. The Cryosphere, 6, 21-34.
# www.the-cryosphere.net/6/21/2012/

# =====================================================================
# Choose where you want the velocity profile to be taken
yp = 15
xp = 15
# =====================================================================

# Script is assuming CISM is using these physical parameter values:
rhoi = 910.0
grav = 9.8101


import sys, os, glob, shutil
from netCDF import *
from math import tan, pi, sin, cos
import numpy as np
import matplotlib.pyplot as plt
from slab import theta, beta, efvs, thickness  # Get the values used to run the experiment

# Calculate scales from Ducowicz paper
eta = beta * thickness / efvs
velscale = rhoi * grav * thickness**2 / efvs


# Get needed variables from the output file
filein = NetCDFFile('slab.out.nc','r')
x1 = filein.variables['x1'][:]
y1 = filein.variables['y1'][:]
level = filein.variables['level'][:]

thk = filein.variables['thk'][:]
if netCDF_module == 'Scientific.IO.NetCDF':
   thk = thk * filein.variables['thk'].scale_factor
topg = filein.variables['topg'][:]
if netCDF_module == 'Scientific.IO.NetCDF':
   topg = topg * filein.variables['topg'].scale_factor
uvel = filein.variables['uvel'][:]
if netCDF_module == 'Scientific.IO.NetCDF':
   uvel = uvel * filein.variables['uvel'].scale_factor


# === Plot the results at the given location ===
# Note we are not plotting like in Fig 3 of paper.
# That figure plotted a profile against zprime.
# It seemed more accurate to plot a profile against z to avoid interpolating model results (analytic solution can be calculated anywhere).

fig = plt.figure(1, facecolor='w')

fig.add_subplot(1,1,1)

# define the scaled x & z variables, with an origin at the bed at this cell
z = level  # elevation / thickness
x = (x1-x1[xp]) / thickness
# calculate rotated zprime coordinates for this column
zprime = -1.0 * x[xp] * sin(theta * pi/180.0) + z * cos(theta * pi/180.0)

# Plot analytic solution for x-component of velocity (eq. 39 in paper)  
plt.plot( sin(theta * pi/180.0) * cos(theta * pi/180.0) * (0.5 * zprime**2 - zprime - 1.0/eta) , \
         z, '-kx', label='Analytic Stokes')
# Plot model results
plt.plot(uvel[0,:,yp,xp] / velscale, z, '--ro', label='CISM') 

plt.legend()
plt.xlabel('Nondimensional velocity')
plt.ylabel('Nondimenstional, unrotated vertical coordinate')
plt.title('Velocity profile at x=' + str(x1[xp]) + ' m, y=' + str(y1[yp]) + ' m')

# Optional plot for comparing analytic solution to Fig. 3 in the paper (model output is not plotted)
#fig = plt.figure(2, facecolor='w')
#plt.plot(level, sin(theta * pi/180.0) * cos(theta * pi/180.0) * (0.5 * level**2 - level - 1.0/eta))
#plt.xlabel("z'")
#plt.ylabel('u')

plt.draw()
plt.show()

filein.close()
