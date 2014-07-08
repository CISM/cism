#!/usr/bin/env python
# This script runs sets up a simple forcing time-series for the dome experiment
# The purpose is for development and testing of time-dependent forcing.
# Currently, it just creates a different acab and artm value for each of 10 years.
# The acab and artm fields are spatially-uniform and negative. 
# It creates the dome.forcing.nc file that is required by the dome.forcing.config file.

# =============================
# Number of time levels to generate
nt=10
# =============================



# Parse options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-c", "--config", dest="configfile", type='string', default='dome.forcincg.config', help="Name of .config file to use to setup and run the dome test case", metavar="FILE")
options, args = optparser.parse_args()

import sys, os, numpy
from netCDF import *
from math import sqrt
from ConfigParser import ConfigParser




# =====================================
# Create a netCDF file according to the information in the config file.
try:
    parser = ConfigParser()
    parser.read(options.configfile)
    nx = int(parser.get('grid','ewn'))
    ny = int(parser.get('grid','nsn'))
    nz = int(parser.get('grid','upn'))
    dx = float(parser.get('grid','dew'))
    dy = float(parser.get('grid','dns'))
    filename = 'dome.forcing.nc'    # This is hard-coded
except:
    sys.exit('Error parsing ' + options.configfile)

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',nt)  # <-- Here is where the number of time levels is set
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = numpy.arange(nt).astype('float32')
netCDFfile.createVariable('x1','f',('x1',))[:] = x
netCDFfile.createVariable('y1','f',('y1',))[:] = y
netCDFfile.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
netCDFfile.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

# Calculate values for the required variables.
#thk  = numpy.zeros([1,ny,nx],dtype='float32')
#topg = numpy.zeros([1,ny,nx],dtype='float32')
#beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')
artm = numpy.zeros([nt,ny,nx],dtype='float32')  # <-- Note the use of nt on these lines
acab = numpy.zeros([nt,ny,nx],dtype='float32')
uvel = numpy.zeros([nt,nz,ny-1,nx-1],dtype='float32')
vvel = numpy.zeros([nt,nz,ny-1,nx-1],dtype='float32')
kinbcmask = numpy.zeros([nt,ny-1,nx-1],dtype='int32')


# Calculate the thickness of the (ellipsoidal) dome of ice
#for i in range(nx):
#  x = float(i-nx/2)/nx
#  for j in range(ny):
#    y = float(j-ny/2)/ny
#    r_squared = x*x+y*y
#    if r_squared < 0.125:
#      thk[0,j,i] = 2000.0 * sqrt(0.125 - r_squared)

# Here is where time-varying values are set
for t in range(nt):
  for j in range(ny):
    acab[t,j,:] = -1.0 * t - j
  artm[t,:,:] = -15.0 - (1.0 * t)
  kinbcmask[t,:,16+t//2:] = 1
  uvel[t,:,:,16+t//2:] = 1.0 + t
  vvel[t,:,:,16+t//2:] = 2.0 + t

#beta[:] = 100000.0
#beta[0,15,1:13] = 1000.0

## specify a sfc temperature field so that temperature evol. can be calc. if desired
#artm[:] = -15.0

# Create the required variables in the netCDF file.
#netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
#netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
#netCDFfile.createVariable('beta','f',('time','y0','x0'))[:] = beta
netCDFfile.createVariable('artm','f',('time','y1','x1'))[:] = artm 
netCDFfile.createVariable('acab','f',('time','y1','x1'))[:] = acab
netCDFfile.createVariable('uvel','f',('time','level','y0','x0'))[:] = uvel[:]
netCDFfile.createVariable('vvel','f',('time','level','y0','x0'))[:] = vvel[:]
netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kinbcmask[:]

netCDFfile.close()

