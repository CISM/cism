#!/usr/bin/env python
# This script runs an experiment with an ice "dome".
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 13, 2010

# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 relase

# Parse options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-c", "--config", dest="configfile", type='string', default='dome.config', help="Name of .config file to use to setup and run the dome test case", metavar="FILE")
optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel', metavar="NUMPROCS")
optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
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
    filename = parser.get('CF input', 'name')
except:
    sys.exit('Error parsing ' + options.configfile)

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',1)
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('staglevel',nz-1)
netCDFfile.createDimension('stagwbndlevel',nz+1)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x
netCDFfile.createVariable('y1','f',('y1',))[:] = y
netCDFfile.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
netCDFfile.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

# Calculate values for the required variables.
thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg = numpy.zeros([1,ny,nx],dtype='float32')
artm = numpy.zeros([1,ny,nx],dtype='float32')
tempstag = numpy.zeros([1,nz+1,ny,nx],dtype='float32')
beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')

# Calculate the thickness of the (ellipsoidal) dome of ice
for i in range(nx):
  x = float(i-nx/2)/nx
  for j in range(ny):
    y = float(j-ny/2)/ny
    r_squared = x*x+y*y
    if r_squared < 0.125:
      thk[0,j,i] = 2000.0 * sqrt(0.125 - r_squared)

# specify a sfc temperature field so that temperature evol. can be calc. if desired
artm[:] = -15.0

# Calculate tempstag and beta values, if desired.  See lines below to enable them being written to file.
tempstag[:] = -0.0
beta[:] = 1.0e8

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
netCDFfile.createVariable('artm','f',('time','y1','x1'))[:] = artm 

# Optional fields that could be added to the initial condition file.  
# Uncomment these lines (and modify their values above), if you want to include them
#netCDFfile.createVariable('tempstag','f',('time','stagwbndlevel','y1','x1'))[:] = tempstag 
#netCDFfile.createVariable('beta','f',('time','y0','x0'))[:] = beta

netCDFfile.close()


# =====================================
# Run CISM
print 'Running CISM'
print '============\n'
if options.parallel == None:
   # Perform a serial run
   os.system(options.executable + ' ' + options.configfile)
else:
   # Perform a parallel run
   if options.parallel <= 0:
      sys.exit( 'Error: Number of processors specified for parallel run is <=0.' )
   else:
      # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
      if os.system('which openmpirun > /dev/null') == 0:
         mpiexec = 'openmpirun -np ' + str(options.parallel)
      elif os.system('which mpirun > /dev/null') == 0:
         mpiexec = 'mpirun -np ' + str(options.parallel)
      elif os.system('which aprun > /dev/null') == 0:
         mpiexec = 'aprun -n ' + str(options.parallel)
      elif os.system('which mpirun.lsf > /dev/null') == 0:
         # mpirun.lsf does NOT need the number of processors (options.parallel)
         mpiexec = 'mpirun.lsf'
      else:
         sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver dome.config')
      runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
      print 'Executing parallel run with:  ' + runstring + '\n\n'
      os.system(runstring)  # Here is where the parallel run is actually executed!


