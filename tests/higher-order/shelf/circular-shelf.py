#!/usr/bin/env python
# This script runs a "Circular Shelf Experiment".
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 9, 2010

# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 relase

import sys, os, glob, shutil, numpy
from netCDF import *
from math import sqrt, exp
from optparse import OptionParser
from ConfigParser import ConfigParser


# Parse command-line options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-c", "--config", dest="configfile", type='string', default='circular-shelf.config', help="Name of .config file to use for the run", metavar="FILE")
optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel [default: perform a serial run]', metavar="NUMPROCS")
optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
optparser.add_option('-b','--smooth-beta',dest='smooth_beta',action='store_true',help='Use a Gaussian function for beta')
optparser.add_option('-d','--dirichlet-center',dest='dirichlet_center',action='store_true',help='Apply Dirichlet boundary condition at the center')
optparser.add_option('-s','--sloped',dest='sloped',action='store_true',help='Use a conically topped ice thickness')
for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()


# Create a netCDF file according to the information in the config file.
parser = ConfigParser()
parser.read(options.configfile)
nx = int(parser.get('grid','ewn'))
ny = int(parser.get('grid','nsn'))
nz = int(parser.get('grid','upn'))
dx = float(parser.get('grid','dew'))
dy = float(parser.get('grid','dns'))
filename = parser.get('CF input', 'name')

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',1)
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x.tolist()
netCDFfile.createVariable('y1','f',('y1',))[:] = y.tolist()
dxs =  dx/2 + x[:-1] # staggered grid
dys =  dy/2 + y[:-1] # staggered grid
netCDFfile.createVariable('x0','f',('x0',))[:] = dxs.tolist()
netCDFfile.createVariable('y0','f',('y0',))[:] = dys.tolist()

# Calculate values for the required variables.
thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg = numpy.empty([1,ny,nx],dtype='float32')
beta = numpy.empty([1,ny-1,nx-1],dtype='float32')
topg[0,:,:] = -2000
# Domain size
Lx = nx*dx
Ly = ny*dy

for i in range(nx):
  x = float(i)/(nx-1) - 0.5   # -1/2 < x < 1/2 
  xx = x*Lx                   # -L/2 < xx < L/2
  for j in range(ny):
    y = float(j)/(ny-1) - 0.5 # -1/2 < y < 1/2
    r = sqrt(x*x+y*y)     # radial distance from the center
    if r < 0.44:          # Inside a circle we have
      thk[0,j,i] = 1000     # constant ice thickness unless
      if options.sloped:  # command line option specifies
        thk[0,j,i] *= (1-r) # conical top
    if options.smooth_beta:  # command line option
      if i < nx-1 and j < ny-1: # beta is on the staggered grid
        yy = y*Ly               # -L/2 < yy < L/2
        beta[0,j,i] = 1.0 + 1.0e10*exp(-(xx*xx+yy*yy)/5.0e5) # Gaussian

if not options.smooth_beta:  # command line option is NOT present
  beta[0,:,:] =  0                             # beta is 0 almost everywhere
  #beta[0,:,:] =  1                             # beta is 1 almost everywhere
  beta[0,ny/2-1:ny/2+1,nx/2-1:nx/2+1] = 1.0e8 # but large in the center

# Add a single bedrock spike in the domain center, to "ground" shelf for 
# bisicles dycore
topg[0,(ny-1)/2-1:(ny-1)/2+2,(nx-1)/2-1:(nx-1)/2+2] = -880. 

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk.tolist()
netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg.tolist()
netCDFfile.createVariable('beta','f',('time','y0','x0'))[:] = beta.tolist()

if options.dirichlet_center:  # command line option
  uvelbc = netCDFfile.createVariable('uvelbc','f',('time','level','y0','x0'))
  vvelbc = netCDFfile.createVariable('vvelbc','f',('time','level','y0','x0'))
  bc = numpy.empty(1,[ny-1,nx-1],dtype='float32')
# boundary condition is NaN almost everywhere
  bc[0,:,:] = float('NaN')
# boundary condition is 0 in the center
  bc[0,ny/2-1:ny/2+2,nx/2-1:nx/2+2] = 0
  for k in range(nz): # loop over levels
    uvelbc[0,k,:,:] = bc.tolist()
    vvelbc[0,k,:,:] = bc.tolist()

netCDFfile.close()


# =====================================
# Run CISM
print 'Running CISM for the circular-shelf experiment'
print '==============================================\n'
if options.parallel == None:
   # Perform a serial run
   runstring = options.executable + ' ' + options.configfile
   print 'Executing serial run with:  ' + runstring + '\n\n'
   os.system(runstring)
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
         sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver circular-shelf.config')
      runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
      print 'Executing parallel run with:  ' + runstring + '\n\n'
      os.system(runstring)  # Here is where the parallel run is actually executed!


