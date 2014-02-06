#!/usr/bin/python
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
import numpy as np
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


# Get the value of flwa specified in the default_flwa parameter in the config file.
# This is the only way this test case supports specifying flwa.
try: 
   flwa=float( subprocess.check_output( 'grep "^default_flwa" ' + configfile + ' | cut -f 3 -d " "', shell='/bin/bash') )
   print 'Parameter used: ' + configfile + ' has specified a flwa value of ' + str(flwa)
except:
   sys.exit('Error: problem getting default_flwa parameter value from the config file')

# Try to get ice density used by the model
try:
   rhoi = float( subprocess.check_output( 'grep "real(dp),parameter :: rhoi =" ../../../libglimmer/glimmer_physcon.F90 | cut -d " " -f 7 | cut -d "." -f 1', shell='/bin/bash' ) )
   print 'Parameter used: ../../../libglimmer/glimmer_physcon.F90 has specified a rhoi value of ' + str(rhoi)
except:
   print 'Warning: problem getting ice density value from ../../../libglimmer/glimmer_physcon.F90  Assuming 910.0 kg/m^3 as a default value.'
   rhoi = 910.0



# open supplied file and get thickness slice needed
filein = NetCDFFile(options.filename,'r')
x1 = filein.variables['x1'][:]
y1 = filein.variables['y1'][:]
time = filein.variables['time'][:]
thk = filein.variables['thk'][:]
if netCDF_module == 'Scientific.IO.NetCDF':
   thk = thk * filein.variables['thk'].scale_factor

# Call the halfar function
thkHalfar = halfarDome(time[timelev]-time[0], x1, y1, flwa, rhoi)

thkDiff = thk[timelev, :, :] - thkHalfar

# Print some stats about the error
print 'Error statistics for cells modeled to have ice:'
print '* Maximum error is ' + str( thkDiff[ np.where( thk[timelev,:,:] > 0.0) ].max() )
print '* Minimum error is ' + str( thkDiff[ np.where( thk[timelev,:,:] > 0.0) ].min() )
print '* Mean error is ' + str( thkDiff[ np.where( thk[timelev,:,:] > 0.0) ].mean() )
print '* Median error is ' + str( np.median(thkDiff[ np.where( thk[timelev,:,:] > 0.0) ]) )

# Plot the results
fig = plt.figure(1, facecolor='w')
markersize = 30.0

fig.add_subplot(1,3,1)
plt.pcolor(x1,y1,thk[timelev,:,:])
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness (m) \n at time ' + str(time[timelev]) )

fig.add_subplot(1,3,2)
plt.pcolor(x1,y1,thkHalfar)
plt.colorbar()
plt.axis('equal')
plt.title('Analytic thickness (m) \n at time ' + str(time[timelev]) )

fig.add_subplot(1,3,3)
plt.pcolor(x1,y1,thkDiff)
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness - Analytic thickness \n at time ' + str(time[timelev]) ) 

plt.draw()
plt.show()

