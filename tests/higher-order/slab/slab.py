#!/usr/bin/env python
# This script runs an experiment with an ice "slab" on an inclined plane.
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Modified from dome.py by Matt Hoffman, Dec. 16, 2013
# Test case described in sections 5.1-2 of:
# J.K. Dukoqicz, 2012. Reformulating the full-Stokes ice sheet model for a more efficient computational solution. The Cryosphere, 6, 21-34.
# www.the-cryosphere.net/6/21/2012/
# However, the implementation used here is based on a more recent unpublished manuscript by J. Dukowicz

# === Physical parameters to adjust, if desired =======================
n = 1 # flow law parameter - only the n=1 case is currently supported
# (implementing the n=3 case would probably require implementing a new efvs option in CISM)
rhoi = 910.0 # kg/m3
grav = 9.1801 # m^2/s

# === Test case parameters to adjust, if desired =======================
theta = 18  # basal inclination angle (degrees)  unpub. man. uses example with theta=18
thickness = 1000.0  # m  thickness in the rotated coordinate system, not in CISM coordinates
# =====================================================================
efvs = 2336041.42829      # hardcoded in CISM for constant viscosity setting (2336041.42829 Pa yr)
# =====================================================================
eta = 10.0   # unpub. man. uses example with eta=10.0
beta = eta / thickness / efvs**-n / (rhoi * grav * thickness)**(n-1)  # Pa yr m^-1
# Note: Fig. 3 in Ducowicz (2013) uses eta=18, where eta=beta*H/efvs
# =====================================================================

if __name__ == '__main__':
  import sys, os, glob, shutil, numpy
  from netCDF import *
  from math import tan, pi, cos
  from ConfigParser import ConfigParser

  # Parse command-line options
  from optparse import OptionParser
  optparser = OptionParser()
  optparser.add_option("-c", "--config", dest="configfile", type='string', default='slab.config', help="Name of .config file to use for the run", metavar="FILE")
  optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel [default: perform a serial run]', metavar="NUMPROCS")
  optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
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
  netCDFfile.createVariable('x1','f',('x1',))[:] = x
  netCDFfile.createVariable('y1','f',('y1',))[:] = y
  netCDFfile.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
  netCDFfile.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

  # Calculate values for the required variables.
  thk  = numpy.zeros([1,ny,nx],dtype='float32')
  topg = numpy.zeros([1,ny,nx],dtype='float32')
  unstagbeta = numpy.zeros([1,ny,nx],dtype='float32')

  # Calculate the geometry of the slab of ice
  thk[:] = thickness / cos(theta * pi/180.0)
  baseElevation = 1000.0 # arbitrary height to keep us well away from sea level
  xmax = x[:].max()
  for i in range(nx):
    topg[0,:,i] = (xmax - x[i]) * tan(theta * pi/180.0) + baseElevation
  offset = 1.0 * float(nx)*dx * tan(theta * pi/180.0)
  parser.set('parameters', 'periodic_offset_ew', str(offset))
  #     Write the new configuration file
  configFile = open(options.configfile,'w')
  parser.write(configFile)
  configFile.close()

  unstagbeta[:] = beta

  # Create the required variables in the netCDF file.
  netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
  netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
  netCDFfile.createVariable('unstagbeta','f',('time','y1','x1'))[:] = unstagbeta

  netCDFfile.close()

  # =====================================
  # Run CISM
  print 'Running CISM for the confined-shelf experiment'
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
           sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver slab.config')
        runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
        print 'Executing parallel run with:  ' + runstring + '\n\n'
        os.system(runstring)  # Here is where the parallel run is actually executed!


