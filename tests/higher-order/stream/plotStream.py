#!/usr/bin/env python
# Script to plot results of stream test case

# Import modules
import sys, os, numpy
from netCDF import *
from runStream import *  # Get all the parameter values and functions defined in the run script
from matplotlib import pyplot as plt

# Parse options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-f", "--file", dest="filename", default="stream.out.nc", help="file to test", metavar="FILE")
for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()



################################
# Start actual script here.

# Open file and get needed vars
f = NetCDFFile(options.filename, 'r')
y0 = f.variables['y0'][:]
ny = y0.shape[0] + 1  # the actual ny dimension, not the size of y0
dy = y0[1]-y0[0]

uvel = f.variables['uvel'][0,:,:,:]  # get all levels at time 0
vvel = f.variables['vvel'][0,:,:,:]  # get all levels at time 0
mintauf = f.variables['tauf'][0,:,:]
#btractx = f.variables['btractx'][0,:,:]
#btracty = f.variables['btracty'][0,:,:]
if netCDF_module == 'Scientific.IO.NetCDF':
   uvel = uvel * f.variables['uvel'].scale_factor
   vvel = vvel * f.variables['vvel'].scale_factor
   mintauf = mintauf * f.variables['tauf'].scale_factor
#   btractx = btractx * f.variables['btractx'].scale_factor
#   btracty = btracty * f.variables['btracty'].scale_factor
#btract = (btractx**2 + btracty**2)**0.5

x0 = f.variables['x0'][:]
xpos = x0.shape[0]/2   # integer division on x-length to get the middle column of the domain

ypos = y0.shape[0]/2   # integer division on y-length to get the middle row of the domain

# Calculate analytic velocity profile - the analytic functions are in runStream.py
if analytic_solution == 'raymond':
    uvel_analytic_profile = raymond_uvel(y0)
    analytic_name = 'raymond analytic solution'
elif analytic_solution == 'schoof':
    uvel_analytic_profile = schoof_uvel(y0)
    analytic_name = 'schoof analytic solution'
else:
    sys.exit("Error: Invalid value for 'analytic_solution'.")

# ===================
# Plot other diagnostics 
# SP: moved this to Fig. 1 so that it lies on bottom)
fig = plt.figure(1, facecolor='w', figsize=(16, 10), dpi=80)

fig.add_subplot(2,2,1)
colors = plt.cm.get_cmap('jet',len(y0))
for j in range(len(y0)):
  plt.plot(x0/1000.0, uvel[ 0,j,:] - uvel[ 0,j,:].mean(), '.-', color=colors(j), label=str(j))
plt.xlabel('distance along flow (km)')
plt.ylabel('surface ALONG flow velocity, demeaned (m/a)')
plt.title('Longit. x-sect. of model vel. for all rows \n(should be constant-valued)')
# kludgy code to get a colorbar for the lines
sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.normalize(vmin=0, vmax=len(y0)))
sm._A = []; 
cb=plt.colorbar(sm)
cb.set_label('y-index')

fig.add_subplot(2,2,2)
for j in range(len(y0)):
  plt.plot(x0/1000.0, uvel[ -1,j,:] - uvel[ -1,j,:].mean(), '.-', color=colors(j), label=str(j))
plt.xlabel('distance along flow (km)')
plt.ylabel('basal ALONG flow velocity, demeaned (m/a)')
#plt.title('Longit. x-sect. of model vel. for all rows (should be constant-valued)')
# kludgy code to get a colorbar for the lines
sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.normalize(vmin=0, vmax=len(y0)))
sm._A = []; 
cb=plt.colorbar(sm)
cb.set_label('y-index')

fig.add_subplot(2,2,3)
for j in range(len(y0)):
  plt.plot(x0/1000.0, vvel[ 0,j,:] - vvel[ 0,j,:].mean(), '.-', color=colors(j), label=str(j))
plt.xlabel('distance along flow (km)')
plt.ylabel('surface ACROSS flow velocity, demeaned (m/a)')
#plt.title('Longit. x-sect. of model vel. for all rows (should be constant-valued)')
# kludgy code to get a colorbar for the lines
sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.normalize(vmin=0, vmax=len(y0)))
sm._A = []
cb=plt.colorbar(sm)
cb.set_label('y-index')

fig.add_subplot(2,2,4)
for j in range(len(y0)):
  plt.plot(x0/1000.0, vvel[-1,j,:] - vvel[-1,j,:].mean(), '.-', color=colors(j), label=str(j))
plt.xlabel('distance along flow (km)')
plt.ylabel('basal ACROSS flow velocity, demeaned (m/a)')
#plt.title('Longit. x-sect. of model vel. for all rows (should be constant-valued)')
# kludgy code to get a colorbar for the lines
sm = plt.cm.ScalarMappable(cmap=colors, norm=plt.normalize(vmin=0, vmax=len(y0)))
sm._A = []
cb=plt.colorbar(sm)
cb.set_label('y-index')

# ===================
# Setup plot of uvel cross-section 
# SP: moved this to Fig. 2 so that it lies on top)
fig = plt.figure(2, facecolor='w', figsize=(12, 10), dpi=80)

fig.add_subplot(2,1,1)
plt.plot(y0/1000.0, uvel_analytic_profile, '-or', label=analytic_name)
plt.plot(y0/1000.0, uvel[ 0,:,xpos], '-xk', label='CISM surface')
plt.plot(y0/1000.0, uvel[-1,:,xpos], '-^k', label='CISM basal')

plt.xlabel('distance across flow (km)')
plt.ylabel('along flow velocity (m/a)')
plt.title(analytic_name +  ' at x=%.1f km'%(x0[xpos]/1000.0))
plt.legend()

fig.add_subplot(2,1,2)
plt.plot(y0/1000.0, mintauf[:, xpos], '-bo', label='Yield stress')
plt.plot(y0/1000.0, numpy.ones(y0.size) * taud, '-k', label='Driving stress')
#plt.plot(y0/1000.0, btract[:, xpos], '-go', label='basal traction')
plt.xlabel('distance across flow (km)')
plt.ylabel('stress (Pa)')
plt.legend()

plt.draw()
plt.show()

#f.close()







