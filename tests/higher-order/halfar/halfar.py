#!/usr/bin/env python
# This script runs an experiment with an ice dome appropriate for the Halfar analytic solution.
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 13, 2010
# Modified for Halfar test case by Matt Hoffman, October 2013.

# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 relase

import sys, os, glob, shutil, numpy
from netCDF import *
from math import sqrt
from ConfigParser import ConfigParser
from halfarDome import halfarDome  # This is located in the current directory
import subprocess

# Check to see if a config file was specified on the command line.
# If not, halfar.config is used.
if len(sys.argv) > 1:
  if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
    print '\nUsage:  python halfar.py [FILE.CONFIG]\n'
    sys.exit(0)
  else:
    configfile = sys.argv[1]
else:
  configfile = 'halfar.config'

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

# Get the value of flwa specified in the default_flwa parameter in the config file.
# This is the only way this test case supports specifying flwa.
try: 
   #flwa=os.system( 'grep "^default_flwa" ' + configfile+ ' | cut -f 3 -d " "')
   flwa=float( subprocess.check_output( 'grep "^default_flwa" ' + configfile+ ' | cut -f 3 -d " "', shell='/bin/bash') )
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


# Calculate the thickness of the halfar dome of ice
thk = halfarDome(0.0, x, y, flwa, rhoi)  # Get the initial time shape from the halfar function
# Note: The halfar solution will assume flwa = 1.0e-16, 
#   so don't modify the default temperature settings.

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg

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
try:
   for files in glob.glob('*.txt')+glob.glob('*.log'):
      # Delete any files already in scratch with these filenames 
      if files in os.listdir('scratch'):
         os.remove(os.path.join('scratch',files))
      # Move the new files to scratch
      shutil.move(files,'scratch')
except:
   print "There was an error in moving the output files to the 'scratch' directory."

print ''
print '================ '
print ''
print 'CISM run completed.  Run ./halfar_results.py to analyze the results.'

