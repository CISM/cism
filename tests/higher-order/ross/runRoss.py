#!/usr/bin/env python
# This script runs a Ross Ice Shelf experiment using Glimmer.
# Output files are written in the "output" subdirectory.
# Two netCDF files are created; one contains "raw" data, the other is an input file for Glimmer.
# After creating the neccessary netCDF input file, Glimmer is run.
# Glimmer writes an additional netCDF output file in the "output" subdirectory.
# When finished, any additional files written by Glimmer are moved to the "scratch" subdirectory.
# After running this script, run plotRoss.py to plot the results.
# See the accompanying README file for more information.
# Written March 18, 2010 by Glen Granzow at the University of Montana.

import os
import glob
import shutil
import numpy
import ConfigParser
from netCDF import *


# =============================================
# === Manually set these options as desired ===
# =============================================
create_files = True
verbose      = False
use_inlets   = (False, True, 'reverse')[1]
# =============================================
# Developer Note: some of these could become command line options below.

# Parse command-line options
from optparse import OptionParser
optparser = OptionParser()
optparser.add_option("-r", "--run", dest="doRun", default=False, action="store_true", help="Including this flag will run CISM.  Excluding it will cause the script to only setup the initial condition file")
optparser.add_option("-c", "--config", dest="configfile", type='string', default='ross.config', help="Name of .config file to use to setup and run the Ross Ice Shelf test case", metavar="FILE")
optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel [default: perform a serial run]', metavar="NUMPROCS")
optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')

for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()


def createArray(data,dtype):
  # nx and ny dimensions are one more than the raw data, 
  # because we are choosing to use the raw data on the velocity grid.
  field = numpy.empty((ny-1,nx-1),dtype=dtype)
  for i in range(nx-1):
    for j in range(ny-1):
      field[j,i] = data[j][i]
  return field

def plot(variable): # used for debugging only
  from matplotlib import pyplot
  pyplot.figure()
  pyplot.imshow(variable,origin='lower',interpolation='nearest')
  pyplot.colorbar()
  pyplot.show()

########## PART I: READ THE INPUT FILES ##########

if create_files:

# Read the main data file into a dictionary mapping names to lists (of lists)
# The dictionary keys are the headers that begin each section in the file
  filename = os.path.join('data','111by147Grid.dat')
  inputfile = open(filename)
  print '\nReading',filename
  currentKey = None
  currentList = list()
  data = dict()
  for line in inputfile:
    line = line.strip()
    if line.startswith('#'):
      if currentKey != None:
        data[currentKey] = currentList
        currentList = list()
      currentKey = line[1:].strip().lower()
    elif len(line) > 0:
      if line.find('.') > 0:
        currentList.append([float(x) for x in line.split()])
      else:
        currentList.append([int(x) for x in line.split()])
  data[currentKey] = currentList
  inputfile.close()

  print 'The',len(data.keys()),'data fields read from 111by147Grid.dat are:\n',data.keys()


  # === Optional ===
  # Manually edit the existency mask so that it extends to be connected to the inlets in inlets.dat 
  data['existency table:'][109][80]=1
  data['existency table:'][97][101]=1
  data['existency table:'][98][100]=1
  data['existency table:'][96][102]=1
  data['existency table:'][96][103]=1
  data['existency table:'][78][121]=1
  data['existency table:'][77][123]=1
  data['existency table:'][79][120]=1
  data['existency table:'][80][120]=1
  data['existency table:'][80][119]=1
  data['existency table:'][51][140]=1
  data['existency table:'][53][138]=1


# Read the kinematic boundary conditions (kbc) mask file.  
# This is a set of (i,j) coordinates specifying where the velocity read from
# the data file should be used as a boundary condition for the model.
  kbc_mask = numpy.zeros((111,147), dtype='i')
  filename = os.path.join('data','kbc.dat')
  inputfile = open(filename)
  print '\nReading',filename
  for line in inputfile:
    i,j = map(int,line.split())
    kbc_mask[i,j] = 1   # 1=where we have Dirichlet; 0=otherwise
  inputfile.close()
  print numpy.sum(kbc_mask),'points were read from kbc.dat'

# Read in the inlets file, which specifies additional Dirichlet boundary conditions
  filename = os.path.join('data','inlets.dat')
  inputfile = open(filename)
  print '\nReading',filename
  counter = [0,0]
  for line in inputfile:
    i, j, azimuth, magnitude = line.split()
    i,j = map(int,(i,j))
    if use_inlets == 'reverse':
      magnitude,azimuth = map(float,(azimuth,magnitude))
    else:
      azimuth,magnitude = map(float,(azimuth,magnitude))
    if use_inlets:
      if verbose:
        indices = '(%d,%d):' % (i,j)
        print 'Changing azimuth  at',indices,data['ice velocity azimuth grid'][i][j],'->',azimuth
        print 'Changing velocity at',indices,data['ice velocity magnitude'][i][j],'->',magnitude
      data['ice velocity azimuth grid'][i][j] = azimuth
      data['ice velocity magnitude'][i][j] = magnitude
    counter[kbc_mask[i,j]] += 1
    print i, j, kbc_mask[i,j] 
    kbc_mask[i,j] += 2
  inputfile.close()
  print 'inlets.dat contains',counter[0],'points that are not in kbc.dat'
  print 'inlets.dat contains',counter[1],'points that are in kbc.dat'

########## PART II: CREATE A NETCDF FILE CONTAINING THE RAW DATA ##########

# This is not necessary, but can be useful to debug the rest of the script, if needed.

#  filename = os.path.join('output','raw.nc')
#  print '\nWriting', filename
#  if netCDF_module == 'netCDF4':
#    netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
#  else:
#    netCDFfile = NetCDFFile(filename,'w')
#  
#  netCDFfile.createDimension('x',data['rows columns number of sub parameters'][0][1])
#  netCDFfile.createDimension('y', data['rows columns number of sub parameters'][0][0])
#  netCDFfile.createVariable('x','f',('x',))[:] = [x[0] for x in data['columns position'][:-1]]
#  netCDFfile.createVariable('y','f',('y',))[:] = [x[0] for x in data['rows position'][:-1]]
#  netCDFfile.createVariable('mask1',    'i',('y','x'))[:] = data['existency table:']
#  netCDFfile.createVariable('azimuth',  'f',('y','x'))[:] = data['ice velocity azimuth grid']
#  netCDFfile.createVariable('velocity', 'f',('y','x'))[:] = data['ice velocity magnitude']
#  netCDFfile.createVariable('thickness','f',('y','x'))[:] = data['thickness']
#  netCDFfile.createVariable('mask2',    'i',('y','x'))[:] = data['reliable velocity obs']
#  netCDFfile.createVariable('seabed',   'f',('y','x'))[:] = data['seabed depth']
#  netCDFfile.createVariable('mask3',    'i',('y','x'))[:] = data['fake ice shelf region']
#  netCDFfile.createVariable('accumulation','f',('y','x'))[:] = data['surface accumulation']
#  netCDFfile.createVariable('bbar',        'f',('y','x'))[:] = data['flowlaw']
#  netCDFfile.createVariable('temperature', 'f',('y','x'))[:] = data['surface temperature']
#  numpy.set_printoptions(threshold='nan')
#  netCDFfile.createVariable('kbc',         'i',('y','x'))[:] = kbc_mask
#  netCDFfile.close()
#  del(netCDFfile) # remove this variable from the name-space (pycdf might fail if we don't)


######### Part II.2  optional plot of kinematic bc positions #######
#  import matplotlib.pyplot as plt
#  plt.imshow(kbc_mask, interpolation='nearest', origin='lower')
#  #plt.imshow(mask1[:,:], interpolation='nearest', origin='lower')
#  for i in range(kbc_mask.shape[1]):
#      for j in range(kbc_mask.shape[0]):
#        if kbc_mask[j,i] == 1:
#          plt.plot(i,j,'og')  # big inlets
#        if kbc_mask[j,i] == 2:
#          plt.plot(i,j,'oc')  # data inlets
#        if data['ice velocity magnitude'][j][i] != 0.0:
#          plt.plot(i,j,'xk')  # nonzero velo
#  plt.colorbar()
#  plt.axis('equal')
#  plt.show()


########## PART III: CREATE THE NETCDF FILE NEEDED BY GLIMMER ##########

# Read the configuration file to get the number of levels to be used
# Also check that ewn, nsn, dew, and dns are correct
  print '\nReading', options.configfile
  configParser = ConfigParser.SafeConfigParser()
  configParser.read(options.configfile)
  nx = int(configParser.get('grid','ewn'))
  ny = int(configParser.get('grid','nsn'))
  nz = int(configParser.get('grid','upn'))
  dx = int(configParser.get('grid','dew'))
  dy = int(configParser.get('grid','dns'))

# nx and ny dimensions are one more than the raw data, 
# because we are choosing to use the raw data on the velocity grid.
  if nx != 147+1:
    print 'WARNING: ewn should be set to 148 in ross.config'
  if ny != 111+1:
    print 'WARNING: nsn should be set to 112 in ross.config'
  if dx != 6822 or dy !=6822:
    print 'WARNING: dew and dns should be set to 6822 in ross.config'

# Put the data into numpy arrays with two extra rows and columns all around
# This reproduces a previous script's output (makerossnc.py)
  mask1        = createArray(data['existency table:'],dtype='i')
  azimuth      = createArray(data['ice velocity azimuth grid'],dtype='f')
  velocity     = createArray(data['ice velocity magnitude'],dtype='f')
  thickness    = createArray(data['thickness'],dtype='f')
#  mask2        = createArray(data['reliable velocity obs'],dtype='i')
  seabed       = createArray(data['seabed depth'],dtype='f')
  mask3        = createArray(data['fake ice shelf region'],dtype='i')
#  accumulation = createArray(data['surface accumulation'],dtype='f')
#  bbar         = createArray(data['flowlaw'],dtype='f')
#  temperature  = createArray(data['surface temperature'],dtype='f')
  kbc          = createArray(kbc_mask,dtype='i')

# Remove any parts of the "fake shelf" that are not in the "existency table"
  if verbose:
    print 'Removing',numpy.sum(numpy.logical_and(mask1==0,mask3!=0)),'points from mask3'
  mask3[mask1==0] = 0

# Set the velocity to zero except where needed as a kinematic boundary condition.
  velocity[kbc == 0] = 0.0
# Get the components of the velocity vector
# Note: velocity1 is vvel, velocity 2 is uvel
  azimuth *= numpy.pi/180.0
  velocity1 = velocity * numpy.cos(azimuth)
  velocity2 = velocity * numpy.sin(azimuth)

# Create the netCDF input file needed by Glimmer
  filename = os.path.join('output','ross.nc')
  print 'Writing', filename
  if netCDF_module == 'netCDF4':
    netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
  else:
    netCDFfile = NetCDFFile(filename,'w')
  netCDFfile.createDimension('time',1)
  netCDFfile.createDimension('x1',nx)
  netCDFfile.createDimension('y1',ny)
  netCDFfile.createDimension('level',nz)
  netCDFfile.createDimension('x0',nx-1) # staggered grid 
  netCDFfile.createDimension('y0',ny-1)
  time  = netCDFfile.createVariable('time','f',('time',))
  x1    = netCDFfile.createVariable('x1','f',('x1',))
  y1    = netCDFfile.createVariable('y1','f',('y1',))
  x0    = netCDFfile.createVariable('x0','f',('x0',))
  y0    = netCDFfile.createVariable('y0','f',('y0',))
  thk   = netCDFfile.createVariable('thk' ,'f',('time','y1','x1'))
  topg  = netCDFfile.createVariable('topg','f',('time','y1','x1'))
  beta  = netCDFfile.createVariable('beta','f',('time','y0','x0'))
  uvel = netCDFfile.createVariable('uvel','f',('time','level','y0','x0'))
  vvel = netCDFfile.createVariable('vvel','f',('time','level','y0','x0'))
  kinbcmask = netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))
  
  time[0] = 0
  x = dx*numpy.arange(nx,dtype='float32')
  y = dx*numpy.arange(ny,dtype='float32')
  x1[:] = x
  y1[:] = y
  x0[:] = dx/2 + x[:-1] # staggered grid
  y0[:] = dy/2 + y[:-1]
  # interpolate thk and topg from the staggered grid onto the scalar grid
  thk[:] = 0.0
  topg[:] = 0.0
  for i in range(1, nx-1):
    for j in range(1, ny-1):
      #assert False
      if numpy.count_nonzero(mask3[j-1:j+1, i-1:i+1]) == 0:
        thk[0,j,i] = thickness[j-1:j+1, i-1:i+1].mean()
      else:
        thk[0,j,i] = 0.0
      topg[0,j,i] = -seabed[j-1:j+1, i-1:i+1].mean()
  # The above line averages the supplied topography from the staggered grid to the unstaggered grid.
  # However, it results in a grounded cell (or maybe cells) in the region specified as the ice shelf.
  # This occurs in a little region just downstream of Roosevelt Island and causes very rapid velocities
  # in that region because it creates a grounded region with a steep slope and beta of 0.
  # Because the test case is only solved for the shelf itself, and there are Dirichlet velocity b.c.
  # at all grounding lines, the topography is not technically needed.  
  # Therefore it is easier to simpy specify a really deep basal topography everywhere to avoid
  # any inadvertent groundings of the ice.
  #topg[:] = -5000.0    
  topg[0,22:42,24:47] -= 500.0

  # extrapolate the edges
  topg[0,0,:] = topg[0,1,:]
  topg[0,-1,:] = topg[0,-2,:]
  topg[0,:,0] = topg[0,:,1]
  topg[0,:,-1] = topg[0,:,-2]
  thk[0,0:2,:] = 0.0  # no ice along bottom two rows of domain to ensure ice shelf front extends to Ross Island side.
  thk[0,-1,:] = thk[0,-2,:]
  thk[0,:,0] = thk[0,:,1]
  thk[0,:,-1] = thk[0,:,-2]
  beta[:] = numpy.zeros((ny-1,nx-1),dtype='float32')
  uvel[:] = numpy.array(nz*[velocity2])
  vvel[:] = numpy.array(nz*[velocity1])
  mask = numpy.logical_and(velocity==0,numpy.logical_or(mask1==1,mask3==1))
  kinbcmask[:] = numpy.int32(numpy.where(mask, 0, 1))

######## Part III.2  optional plot of kinematic bc positions #######
#  import matplotlib.pyplot as plt
#  fig = plt.figure(2, facecolor='w', figsize=(10, 4), dpi=100)
#
#  #plt.imshow(kinbcmask[0,:,:], interpolation='nearest', origin='lower')
#  plt.imshow(mask1[:,:], interpolation='nearest', origin='lower')
#  for i in range(kbc.shape[1]):
#      for j in range(kbc.shape[0]):
#        if kbc[j,i] == 1:
#          plt.plot(i,j,'og')  # big inlets
#        if kbc[j,i] == 2:
#          plt.plot(i,j,'oc')  # data inlets
#        if velocity[j,i] != 0.0:
#          plt.plot(i,j,'xk')  # nonzero velo
#  plt.colorbar()
#  plt.axis('equal')
#  plt.show()

  netCDFfile.close()



########## PART IV: RUN GLIMMER ##########
if options.doRun:
    # =====================================
    # Run CISM
    print 'Running CISM for the Ross Ice Shelf experiment'
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
             sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver ross.config')
          runstring = mpiexec + ' ' + options.executable + ' ' + options.configfile
          print 'Executing parallel run with:  ' + runstring + '\n\n'
          os.system(runstring)  # Here is where the parallel run is actually executed!
    # Clean up by moving .log file for this run to the 'output' directory 
    try:
      shutil.move(options.configfile+'.log', './output')
    except:
      pass  # if there was a problem, don't worry about it...
else:
    print "\nInitial condition file " + filename + " has been created, but CISM has not been run.  Either run the model manually or use 'python runRoss.py --help' for details of how to run the model with this script."


