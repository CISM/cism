#!/usr/bin/env python
# A script to compare CISM model output to the Halfar analytic solution of SIA evolution of a dome.
# Matt Hoffman, LANL, October 2013

import sys
import datetime
import subprocess
try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 '
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')
from optparse import OptionParser
from ConfigParser import ConfigParser
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from halfarDome import halfarDome   # located in current directory

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to test", metavar="FILE")
parser.add_option("-t", "--time", dest="t", help="which time level to use", metavar="T")
parser.add_option("-c", "--config", dest="configfile", help="name of config file used for model run", metavar="CONFIG")

options, args = parser.parse_args()
if not options.filename:
   options.filename = 'halfar.out.nc'
   print 'No file specified.  Attempting to use halfar.out.nc'
if options.t:
   timelev = int(options.t)
else:
   timelev = -1
   print 'No time level specified.  Attempting to use final time.'
if options.configfile:
   configfile = options.configfile
else:
   configfile = 'halfar.config'
   print 'No config file specified.  Attempting to use halfar.config.'

# Open config file for reading
conparser = ConfigParser()
conparser.read(configfile)

# Get the value of flwa specified in the default_flwa parameter in the config file.
# This is the only way this test case supports specifying flwa.
try: 
   flwa = float(conparser.get('parameters','default_flwa'))
   print 'Parameter used: ' + configfile + ' has specified a flwa value of ' + str(flwa)
   flow_law = int(conparser.get('options','flow_law'))
   if flow_law != 0:
      sys.exit('Error: The option "flow_law" must be set to 0 for the test case to work properly.')
except:
   sys.exit('Error: problem getting default_flwa parameter value from the config file')

# Try to get ice density used by the model
try:
   rhoi = float( subprocess.check_output( 'grep "real(dp),parameter :: rhoi =" ../../libglimmer/glimmer_physcon.F90 | cut -d " " -f 7 | cut -d "." -f 1', shell='/bin/bash' ) )
   print 'Parameter used: ../../libglimmer/glimmer_physcon.F90 has specified a rhoi value of ' + str(rhoi)
except:
   print 'Warning: problem getting ice density value from ../../../libglimmer/glimmer_physcon.F90  Assuming 910.0 kg/m^3 as a default value.'
   rhoi = 910.0



# open supplied file and get thickness slice needed
filein = NetCDFFile(options.filename,'r')
x1 = filein.variables['x1'][:]
dx = x1[1]-x1[0]
y1 = filein.variables['y1'][:]
ny = y1.size
time = filein.variables['time'][:]
thk = filein.variables['thk'][:]
if netCDF_module == 'Scientific.IO.NetCDF':
   thk = thk * filein.variables['thk'].scale_factor

# Call the halfar function
thkHalfar = halfarDome(time[timelev]-time[0], x1, y1, flwa, rhoi)

thkDiff = thk[timelev, :, :] - thkHalfar
thkDiffIce = thkDiff[ np.where( thk[timelev,:,:] > 0.0) ]  # Restrict to cells modeled to have ice
RMS = ( (thkDiffIce**2).sum() / float(len(thkDiffIce)) )**0.5

# Print some stats about the error
print '\nError statistics for cells modeled to have ice (in m):'
print '* RMS error = ' + str( RMS )
print '* Maximum error is ' + str( thkDiffIce.max() )
print '* Minimum error is ' + str( thkDiffIce.min() )
print '* Mean error is ' + str( thkDiffIce.mean() )
print '* Median error is ' + str( np.median(thkDiffIce) )
print '* Mean absolute error = ' + str( np.absolute(thkDiffIce).mean() )
print '* Median absolute error = ' + str( np.median(np.absolute(thkDiffIce)) )
print ''


# =========================
# Plot the results
fig = plt.figure(1, facecolor='w', figsize=(10, 4), dpi=100)
gray = np.ones(3)*0.8

fig.add_subplot(1,3,1)
plt.pcolor(x1/1000.0,y1/1000.0,  ma.masked_values(thk[timelev,:,:], 0.0) )
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness (m) \n at time ' + str(time[timelev]) )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

fig.add_subplot(1,3,2)
plt.pcolor(x1/1000.0,y1/1000.0, ma.masked_values(thkHalfar, 0.0) )
plt.colorbar()
plt.axis('equal')
plt.title('Analytic thickness (m) \n at time ' + str(time[timelev]) )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

fig.add_subplot(1,3,3)
plt.pcolor(x1/1000.0,y1/1000.0, ma.masked_values(thkDiff, 0.0))
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness - Analytic thickness \n at time ' + str(time[timelev]) ) 
plt.xlabel('x (km)'); plt.ylabel('y (km)')


# =========================
# optional second figure - cross section through center of dome

#fig = plt.figure(2, facecolor='w', figsize=(10, 4), dpi=100)
#yind = ny//2
#print yind, y1[yind]

#x1dense = np.linspace(x1[0], x1[-1], 1000)
#thkHalfarDense = halfarDome(time[timelev]-time[0], x1dense, y1[[0, yind, -1]], flwa, rhoi)

#plt.step(x1/1000.0 + 0.5*dx/1000.0, thk[timelev,yind,:], '-r', label='model')
#plt.plot(x1/1000.0, thk[timelev,yind,:], '.r')

#plt.plot(x1dense/1000.0, thkHalfarDense[1,:], '-k', label='analytic')

#plt.stem(x1/1000.0, (thk[timelev,yind,:] - thkHalfar[yind,:]) * 10.0, '-.b', label='error x10.0')

#plt.xlabel('x (km)'); plt.ylabel('Elevation (m)')
#plt.legend()



plt.draw()
plt.show()

