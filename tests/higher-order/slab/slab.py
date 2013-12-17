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

# === Test case parameters to adjust, if desired =======================
theta = 18  # basal inclination angle (degrees)
thickness = 1000.0  # m  thickness in the rotated coordinate sytem, not in CISM coordinates
# =====================================================================
efvs = 1.0e7      # hardcoded in CISM (10^7 Pa yr)
# =====================================================================
beta = 18.0 / thickness * efvs  # Pa yr m^-1
# Note: Fig. 3 in Ducowicz (2013) uses eta=18, where eta=beta*H/efvs
# =====================================================================

if __name__ == '__main__':
  import sys, os, glob, shutil, numpy
  from netCDF import *
  from math import tan, pi, cos
  from ConfigParser import ConfigParser


  # Check to see if a config file was specified on the command line.
  # If not, slab.config is used.
  if len(sys.argv) > 1:
    if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
      print '\nUsage:  python slab.py [FILE.CONFIG]\n'
      sys.exit(0)
    else:
      configfile = sys.argv[1]
  else:
    configfile = 'slab.config'

  # Check to see if #procs specified, relevant when running the code in parallel. 
  # If not, serial run (#procs==1) is performed. To run in parallel, the configure
  # file must be specifed, but the nu,ber of processors does not
  if len(sys.argv) > 2:
      nprocs = sys.argv[2]
  else:
    nprocs = '1'

  # Create a netCDF file according to the information in the config file.
  parser = ConfigParser()
  parser.read(configfile)
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
  for i in range(nx):
    topg[0,:,i] = x[i] * tan(theta * pi/180.0) + baseElevation
  offset = -1.0 * float(nx)*dx * tan(theta * pi/180.0)
  parser.set('parameters', 'periodic_offset_ew', str(offset))
  #     Write the new configuration file
  configFile = open(configfile,'w')
  parser.write(configFile)
  configFile.close()

  unstagbeta[:] = beta

  # Create the required variables in the netCDF file.
  netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
  netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
  netCDFfile.createVariable('unstagbeta','f',('time','y1','x1'))[:] = unstagbeta

  netCDFfile.close()

  # Run Glimmer
  print 'Running Glimmer/CISM'
  if len(sys.argv) > 2:
     #os.system('aprun -n'+nprocs+' ./simple_glide '+configfile+'')  # support for MPI runs is here
     os.system('mpirun -np '+nprocs+' ./simple_glide '+configfile+'') 
  else:
     os.system('echo '+configfile+' | ./simple_glide')

  # Clean up by moving extra files written by Glimmer to the "scratch" subdirectory
  # Look for files with extension "txt", "log", or "nc"
  for files in glob.glob('*.txt')+glob.glob('*.log'):
  # Delete any files already in scratch with these filenames 
    if files in os.listdir('scratch'):
      os.remove(os.path.join('scratch',files))
  # Move the new files to scratch
    shutil.move(files,'scratch')
