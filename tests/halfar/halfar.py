#!/usr/bin/env python
# This script runs an experiment with an ice dome appropriate for the Halfar analytic solution.
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Modified from dome.py script written by Glen Granzow at the University of Montana on April 13, 2010
# Modified for Halfar test case by Matt Hoffman, October 2013.

# Parse options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-c", "--config", dest="configfile", type='string', default='halfar.config', help="Name of .config file to use to setup and run the Halfar test case.", metavar="FILE")
optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel.  If absent, the run will be in serial (which is required for Glide SIA dycore).', metavar="NUMPROCS")
optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()
options, args = optparser.parse_args()

import sys, os, numpy
from netCDF import *
from math import sqrt
from ConfigParser import ConfigParser
from halfarDome import halfarDome  # This is located in the current directory
import subprocess


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

print 'Creating', filename
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
   flwa = float(parser.get('parameters','default_flwa'))
   print 'Parameter used: ' + options.configfile + ' has specified a flwa value of ' + str(flwa)
   flow_law = int(parser.get('options','flow_law'))
   if flow_law != 0:
      sys.exit('Error: The option "flow_law" must be set to 0 for the test case to work properly.')
except:
   sys.exit('Error: problem getting default_flwa parameter value from the config file ' + options.configfile)

# Try to get ice density used by the model
try:
   rhoi = float( subprocess.check_output( 'grep "real(dp),parameter :: rhoi =" ../../libglimmer/glimmer_physcon.F90 | cut -d " " -f 7 | cut -d "." -f 1', shell='/bin/bash' ) )
   print 'Parameter used: ../../libglimmer/glimmer_physcon.F90 has specified a rhoi value of ' + str(rhoi)
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
         sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./simple_glide dome.config')
      runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
      print 'Executing parallel run with:  ' + runstring + '\n\n'
      os.system(runstring)  # Here is where the parallel run is actually executed!


print ''
print '================ '
print ''
print 'CISM run completed.  Run "python halfar_results.py" to analyze the results.'

