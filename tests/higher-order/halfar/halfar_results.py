#!/usr/bin/python
# A script to compare CISM model output to the Halfar analytic solution of SIA evolution of a dome.
# Matt Hoffman, LANL, October 2013

import sys
import datetime
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

options, args = parser.parse_args()
if not options.filename:
   options.filename = 'halfar.out.nc'
   print 'No file specified.  Attempting to use halfar.out.nc'
if options.t:
   timelev = int(options.t)
else:
   timelev = -1
   print 'No time level specified.  Attempting to use final time.'



# open supplied file and get thickness slice needed
filein = NetCDFFile(options.filename,'r')
x1 = filein.variables['x1'][:]
y1 = filein.variables['y1'][:]
time = filein.variables['time'][:]
thk = filein.variables['thk'][:]
if netCDF_module == 'Scientific.IO.NetCDF':
   thk = thk * filein.variables['thk'].scale_factor

# Call the halfar function
thkHalfar = halfarDome(time[timelev]-time[0], x1, y1)

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

