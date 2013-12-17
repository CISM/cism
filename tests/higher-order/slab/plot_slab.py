#!/usr/bin/env python
# This script plots result of an experiment with an ice "slab" on an inclined plane.
# Written by Matt Hoffman, Dec. 16, 2013
# Test case described in sections 5.1-2 of:
# J.K. Dukoqicz, 2012. Reformulating the full-Stokes ice sheet model for a more efficient computational solution. The Cryosphere, 6, 21-34.
# www.the-cryosphere.net/6/21/2012/



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
x0 = filein.variables['x0'][:]
y0 = filein.variables['y0'][:]
level = filein.variables['level'][:]

# =====================================================================
# Choose where you want the velocity profile to be taken
# use integer floor division operator to get an index close to the center 
xp = len(x0)//2
yp = len(y0)//2
#yp = 15
#xp = 15
# =====================================================================
print 'Using x index of '+str(xp)+'='+str(x0[xp])
print 'Using y index of '+str(yp)+'='+str(y0[yp])

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

fig = plt.figure(1, facecolor='w', figsize=(12, 6))

# define the scaled x & z variables, with an origin at the bed at this cell
z = level  # elevation / thickness
x = (x0-x0[xp]) / thickness
# calculate rotated zprime coordinates for this column
zprime = -1.0 * x[xp] * sin(theta * pi/180.0) + (1.0-z) * cos(theta * pi/180.0)
# Calculate analytic solution for x-component of velocity (eq. 39 in paper) for the CISM-column
uvelAnalyticScaled =  sin(theta * pi/180.0) * cos(theta * pi/180.0) * (0.5 * zprime**2 - zprime - 1.0/eta)

### 1. Plot as nondimensional variables
# Plot analytic solution 
fig.add_subplot(1,2,1)
plt.plot(uvelAnalyticScaled, z, '-kx', label='Analytic Stokes')
# Plot model results
plt.plot(uvel[0,:,yp,xp] / velscale, z, '--ro', label='CISM') 
plt.gca().invert_yaxis()  # put 0.0 (surface) on top & 1.0 (bed) on bottom
plt.legend()
plt.xlabel('Nondimensional velocity')
plt.ylabel('Nondimenstional, unrotated vertical coordinate')
plt.title('Velocity profile at x=' + str(x0[xp]) + ' m, y=' + str(y0[yp]) + ' m\n(Scaled coordinates)')

### 2. Plot as dimensional variables
# Plot analytic solution for x-component of velocity (eq. 39 in paper)  
fig.add_subplot(1,2,2)
plt.plot(uvelAnalyticScaled * velscale, (1.0-z) * thk[0,yp,xp] + topg[0,yp,xp], '-kx', label='Analytic Stokes')
# Plot model results
plt.plot(uvel[0,:,yp,xp], (1.0-z) * thk[0,yp,xp] + topg[0,yp,xp], '--ro', label='CISM') 
plt.legend()
plt.xlabel('velocity (m/yr)')
plt.ylabel('elevation (m)')
plt.title('Velocity profile at x=' + str(x0[xp]) + ' m, y=' + str(y0[yp]) + ' m\n(Unscaled coordinates)')

#################
# Now plot maps to show if the velocities vary over the domain (they should not)
fig = plt.figure(2, facecolor='w', figsize=(12, 6))
fig.add_subplot(1,2,1)
uvelDiff = uvel[0,0,:,:] - uvel[0,0,yp,xp]
tol = 1.0e-12
#plt.pcolor(x0,y0,uvelDiff,vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol)
plt.imshow(uvelDiff,extent=(x0.min(), x0.max(), y0.min(), y0.max()),vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol, interpolation="none")
plt.xlabel('x position (m)')
plt.ylabel('y position (m)')
plt.colorbar()
plt.plot(x0[xp],y0[yp],'k*',markersize=9)
plt.title('Map of difference of surface x-velocity from\nvalue at reference location (m/yr)')
plt.axis('equal')

fig.add_subplot(1,2,2)
uvelDiff = uvel[0,-1,:,:] - uvel[0,-1,yp,xp]
tol = 1.0e-12
#plt.pcolor(x0,y0,uvelDiff,vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol)
plt.imshow(uvelDiff,extent=(x0.min(), x0.max(), y0.min(), y0.max()),interpolation="none",vmin=uvelDiff.min()-tol, vmax=uvelDiff.max()+tol)
plt.xlabel('x position (m)')
plt.ylabel('y position (m)')
plt.colorbar()
plt.plot(x0[xp],y0[yp],'k*',markersize=9)
plt.title('Map of difference of basal x-velocity from\nvalue at reference location (m/yr)')
plt.axis('equal')


# Optional plot for comparing analytic solution to Fig. 3 in the paper (model output is not plotted)
#fig = plt.figure(9, facecolor='w')
#plt.plot(level, sin(theta * pi/180.0) * cos(theta * pi/180.0) * (0.5 * level**2 - level - 1.0/eta))
#plt.xlabel("z'")
#plt.ylabel('u')

plt.draw()
plt.show()

filein.close()
