#!/usr/bin/env python

# This script sets up initial conditions for the MISMIP3d experiment.
# See this paper for details about the MISMIP3d experiment:
# Pattyn et al., Grounding-line migration in plan-view marine ice-sheet models: results of the ice2sea MISMIP3d intercomparison, J. Glaciol., 59, doi:10.3189/2013JoG12J129, 2013.

import sys, os
import shutil
import fileinput
import numpy as np
from netCDF4 import Dataset
from ConfigParser import ConfigParser
from optparse import OptionParser


#############
# Constants #
#############


xDomain = 800000.0     # domain x-dimension (m)
yDomain = 100000.0     # domain y-dimension (m)
initThickness = 500.   # initial uniform ice thickness
accum = 0.5            # uniform accumulation (m/yr)


####################################
# Function used later in the code #
####################################


# The following function returns the linear bed topography as in Pattyn et al. (2013).
def computeBed(x):
    # Input x has to be in km.
    slope = -1.
    b0    = -100.    # m
    eps_b = 1e-10    # small regularization number
    abs_x = np.sqrt(x**2 + eps_b**2)  # regularizing to avoid problems at the divide
  
    b = b0 + slope*abs_x
    return b



########
# Code #
########


# Parse options
optparser = OptionParser()

optparser.add_option('-c', '--config', dest='configfile',   type='string',  default='mismip3d.config.template', help='config file name for setting up the MISMIP3d experiment', metavar='FILE')
optparser.add_option('-e', '--exec',   dest='executable',                   default='cism_driver',help='Set path to the CISM executable', metavar='EXECUTABLE')
optparser.add_option('-x', '--expt',   dest='experiment',    type='string', default = 'all',    help='MISMIP3d experiment to set up', metavar='EXPT')
optparser.add_option('-t', '--tstep',  dest='timestep',      type='float',  default = 1,        help='time step (yr)', metavar='TSTEP')
optparser.add_option('-r', '--res',    dest='resolution',    type='int',    default = 2000,     help='grid resolution (m)', metavar='RES')
optparser.add_option('-v', '--vlevel', dest='vertlevels',    type='int',    default = 3,        help='no. of vertical levels', metavar='VLEVEL')
optparser.add_option('-a', '--approx', dest='approximation', type='string', default = 'DIVA',   help='Stokes approximation (SSA, DIVA, BP)', metavar='APPROXIMATION')
optparser.add_option('-b', '--basal',  dest='basalFriction', type='string', default='powerlaw', help='Basal friction law (powerlaw, schoof)', metavar='BASALFRICTION')
optparser.add_option('-y', '--year',   dest='yearsStnd',     type='int',    default = 20000,    help='Length of Stnd run (yr)', metavar='YEARSPINUP')

optparser.add_option

for option in optparser.option_list:
    if option.default != ('NO', 'DEFAULT'):
        option.help += (' ' if option.help else '') + '[default: %default]'
options, args = optparser.parse_args()

if options.experiment == 'all':
    experiments = ['Stnd','P75S','P75R']
    print 'Setting up all the MISMIP3d experiments'
elif options.experiment == 'allP75':
    experiments = ['P75S','P75R']
    print 'Setting up P75S and P75R experiments'
elif options.experiment in ['Stnd','P75S','P75R']:
    experiments = [options.experiment]
    print 'Setting up experiment', options.experiment
else:
    sys.exit('Please specify experiment(s) from this list: Stnd, P75S, P75R')


# If there is not already a link to cism_driver in the main directory, then make one.
# Each subdirectory will link to cism_driver in the main directory.
if options.executable != 'cism_driver':
    # Remove the existing link, if present.
    os.unlink('cism_driver')
    # Make the new link.
    os.symlink(options.executable, 'cism_driver')


# Set grid resolution.
if options.resolution == 8000:
    dx = 8000.0
    dy = 8000.0
elif options.resolution == 4000:
    dx = 4000.0
    dy = 4000.0
elif options.resolution == 2000:
    dx = 2000.0
    dy = 2000.0
elif options.resolution == 1000:
    dx = 1000.0
    dy = 1000.0
elif options.resolution == 500:
    dx = 500.0
    dy = 500.0
elif options.resolution == 250:
    dx = 250.0
    dy = 250.0
else:
    sys.exit('Please choose from among the following resolutions (m): 8000, 4000, 2000, 1000, 500, 250')

if options.vertlevels >= 2:
    nz = options.vertlevels
else:
    sys.exit('Error: must have at least 2 vertical levels')

print 'MISMIP3d grid resolution (m) =', options.resolution
print 'Number of vertical levels =', nz

# Set number of grid cells in each direction.
# Include a few extra cells in the x direction to handle boundary conditions.
nx = int(xDomain/dx) + 4
ny = int(yDomain/dy)

# Copy the config template to a new master config file.
masterConfigFile = 'mismip3d.config'

try:
    shutil.copy(options.configfile, masterConfigFile)
except OSError:
    sys.exit('Could not copy', options.configfile)

print 'Creating master config file', masterConfigFile

# Read the master config file.
config = ConfigParser()
config.read(masterConfigFile)

# Set the grid variables in the master config file.
config.set('grid', 'ewn', nx)
config.set('grid', 'nsn', ny)
config.set('grid', 'upn', nz)
config.set('grid', 'dew', dx)
config.set('grid', 'dns', dy)

# Set the time step in the msster config file.
# Set the diagnostic interval to the same value (not necessary, but helpful for debugging).

config.set('time', 'dt', options.timestep)
config.set('time', 'dt_diag', options.timestep)

# Set Stokes approximation in config file.
if options.approximation == 'SSA':
    which_ho_approx = 1
    print 'Using SSA velocity solver'
elif options.approximation == 'DIVA':
    which_ho_approx = 4
    print 'Using DIVA velocity solver'
elif options.approximation == 'BP':
    which_ho_approx = 2
    print 'Using Blatter-Pattyn velocity solver'
else:
    which_ho_approx = 4
    print 'Defaulting to DIVA velocity solver'

config.set('ho_options', 'which_ho_approx', which_ho_approx)

# Config settings related to basal friction law.
# Note: Each of these friction laws is associate with certain basal parameters.
#       The desired parameters should be set in the config template.
if options.basalFriction == 'powerlaw':
    p_ocean_penetration = 0
    print 'Using basal friction power law (Schoof2007)'
elif options.basalFriction == 'schoof':
    p_ocean_penetration = 1
    print 'Using Schoof 2005 basal friction law'
else:
    p_ocean_penetration = 0   #  is default
    print 'Defaulting to Powerlaw basal friction law'

config.set('parameters', 'p_ocean_penetration', p_ocean_penetration)

# Config setting related to spin up time.
yearsStnd = float(options.yearsStnd)
config.set('time', 'tend', yearsStnd)

# Write to the master config file.
with open(masterConfigFile, 'w') as configfile:
    config.write(configfile)


print 'years of Stnd experiment =', yearsStnd
restartfreqStnd = min(1000.0, options.yearsStnd)    # can be changed by the user if needed
print 'Stnd restart frequency =', restartfreqStnd


# Create the netCDF input file according to the information in the config file.
try:
    parser = ConfigParser()
    parser.read(options.configfile)
    initfile = parser.get('CF input', 'name')
except OSError:
    sys.exit('Error parsing ' + options.configfile)
    
print 'Creating input file', initfile
ncfile = Dataset(initfile, 'w')
    

# Create dimensions.
# Note: (x0,y0) = staggered (velocity) grid.
#       (x1,y1) = unstaggered (scalar) grid.
ncfile.createDimension('time',1)
ncfile.createDimension('x1',nx)
ncfile.createDimension('y1',ny)
ncfile.createDimension('x0',nx-1)
ncfile.createDimension('y0',ny-1)
ncfile.createDimension('level',nz)
ncfile.createDimension('staglevel',nz-1)
ncfile.createDimension('stagwbndlevel',nz+1)  # similar to staglevel but including boundaries

# Create time and grid variables.
# Note: (x1,y1) are loadable and need to be in the input file.
#       (x0,y0) are not loadable, but are derived in CISM from (x1,y1). May not be needed.
ncfile.createVariable('time','f4',('time',))[:] = [0]
x1 = ncfile.createVariable('x1', 'f4', ('x1',))
y1 = ncfile.createVariable('y1', 'f4', ('y1',))
x0 = ncfile.createVariable('x0', 'f4', ('x0',))
y0 = ncfile.createVariable('y0', 'f4', ('y0',))

# Create 2D input fields.
thk  = ncfile.createVariable('thk',  'f4', ('time','y1','x1'))
topg = ncfile.createVariable('topg', 'f4', ('time','y1','x1'))
acab = ncfile.createVariable('acab', 'f4', ('time','y1','x1'))
uvel = ncfile.createVariable('uvel', 'f4', ('time','level','y0','x0'))
vvel = ncfile.createVariable('vvel', 'f4', ('time','level','y0','x0'))
C_space_factor = ncfile.createVariable('C_space_factor','f4',('time','y1','x1'))
kinbcmask      = ncfile.createVariable('kinbcmask',     'i4', ('time','y0','x0'))  # kinematic BC mask

# Compute x and y on each grid.
# Note: (1) The x origin is placed at the center of the second cell from the left.
#           This assumes that kinbcmask = 1 at the first vertex from the left.
#           Thus the left edge of the grid has x = -3*dx/2.
#       (2) The y origin is placed at the bottom edge of the CISM grid.
#           The line of central symmetry runs along cell edges at y = 40 km.

x = dx*np.arange(nx,dtype='float32')  # x = 0, dx, 2*dx, etc.
y = dy*np.arange(ny,dtype='float32')  # y = 0, dy, 2*dy, etc.

x1[:] = x[:] - dx         # x1 = -dx, 0, dx, ..., (nx-2)*dx - dx/2
y1[:] = y[:] + dy/2.      # y1 = dy/2, 3*dy/2, ..., (ny-1)*dy - dy/2

x0[:] = x[:-1] - dx/2.    # x0 = -dx/2, dx/2, 3*dx/2, ..., (nx-2)*dx
y0[:] = y[:-1] + dy       # y0 = dy, 2*dy, ..., (ny-1)*dy


# Set bed topography.
for i in range(nx):
    topg[:,:,i] = computeBed(x1[i]/1.e3) # x1 is in [m] and we input in [km]

# Set initial thickness.
thk[0,:,:] = initThickness

# Set the surface mass balance.
acab[:,:,:] = accum

# Set initial velocity to zero (probably not necessary).
uvel[:,:,:,:] = 0.
vvel[:,:,:,:] = 0.

# Set kinematic velocity mask.
# Where kinbcmask = 1, the velocity is fixed at its initial value.
# Note: Although there is no ice on the RHS of the domain, we need kinbcmask =1 there
#       to preserve symmetry with the LHS (since east-west BCs are formally periodic).
kinbcmask[:,:,:]  = 0   # initialize to 0 everywhere
kinbcmask[:,:,0]  = 1   # mask out left-most column
kinbcmask[:,:,-1] = 1   # mask out right-most column

# Set the spatially variable basal friction coefficient.
C_space_factor[0,:,:] = 1.

ncfile.close()

print 'Experiments:', experiments

# Loop through experiments.
for expt in experiments:
    
    # For each experiment, make a suitable config file and set up a subdirectory.
    print 'Creating config file for experiment', expt
    
    # Make a copy of the mismip3dInit config file.
    # Below, this copy will be tailored for the chosen MISMIP3d experiment,
    # without changing the settings used for the Stnd experiment.
    
    newConfigFile = 'mismip3d' + expt + '.config'
    print 'Config file for this experiment:', newConfigFile
    shutil.copy(masterConfigFile, newConfigFile)
    
    # Read the new config file.
    config = ConfigParser()
    config.read(newConfigFile)

# Experiment-specific settings.

    if expt == 'Stnd':
        tstart      = 0.0
        tend        = yearsStnd
        inputdir    = '../'
        inputfile   = initfile
        inputslice  = 1
        outputfreq  = min(1000.0, restartfreqStnd)
        restartfreq = restartfreqStnd
    elif expt == 'P75S':
        tstart       = 0.0
        tend         = 100.0
        inputdir     = '../Stnd/'
        inputfile    = 'mismip3dStnd.restart.nc'
        inputfileOut = 'mismip3dStnd.out.nc'  # we will need f_ground from this file to calculate the GL (not in the restart file)
        inputslice   = int(yearsStnd/restartfreqStnd)
        outputfreq   = 10.0
        restartfreq  = 100.0
    elif expt == 'P75R':
        tstart      = 0.0
        tend        = yearsStnd
        inputdir    = '../P75S/'
        inputfile   = 'mismip3dP75S.restart.nc'
        inputslice  = 1
        outputfreq  = 1000.0
        restartfreq = 1000.0

    # Set the start and end times
    config.set('time', 'tstart', tstart)
    config.set('time', 'tend',   tend)
    
    # Change the default comment
    comment = 'MISMIP3d experiment ' + expt
    config.set('CF default', 'comment', comment)
    
    # Set input file and time slice.
    # Note: This method may not be robust for Stnd and P75R runs that start and restart.
    #       For this reason, the script mismip3dRun.py makes sure the 'time' entry
    #       in [CF input] corresponds to the final time slice.
    print 'Input file:', inputfile
    if expt=='P75S':
       config.set('CF input', 'nameOut', inputfileOut)

    config.set('CF input', 'name', inputfile)
    config.set('CF input', 'time', inputslice)

    # Set the output filename in the section [CF output1].
    outfilename = 'mismip3d' + expt + '.out.nc'
    print 'Output file:', outfilename
    config.set('CF output', 'name',      outfilename)
    config.set('CF output', 'frequency', outputfreq)

    # Set restart info.  This should be in the section called '[CF output]'.
    # Note: Each experiment (except Stnd) writes out only one time slice to a restart file.
    restartfilename = 'mismip3d' + expt + '.restart.nc'
    print 'Restart file:', restartfilename
    config.set('CF restart', 'name',      restartfilename)
    config.set('CF restart', 'variables', 'restart')
    config.set('CF restart', 'frequency', restartfreq)
    config.set('CF restart', 'start',     tstart + restartfreq)
    config.set('CF restart', 'xtype',     'double')

    # Write to the new config file.
    with open(newConfigFile, 'w') as configfile:
        config.write(configfile)

    # Create a subdirectory named for the experiment, and stage the run there.
    try:
        os.mkdir(expt)
        print 'Created subdirectory', expt
    except OSError:
        print 'Subdirectory', expt, 'already exists'

    os.chdir(expt)
    
    # Move the config file from the parent directory to the subdirectory.
    shutil.move('../' + newConfigFile, newConfigFile)
    print 'Created config file', newConfigFile
    
    # Link to the cism_driver executable in the parent directory.
    try:
        os.symlink('../cism_driver', 'cism_driver')
    except OSError:
        pass   # link to cism_driver already exists

    # Link to the input file in the appropriate directory.
    try:
        os.symlink(inputdir + inputfile, inputfile)
    except OSError:
        pass   # link to the input file already exists

    if expt=='P75S':
        try:
           os.symlink(inputdir + inputfileOut, inputfileOut)
        except OSError:
           pass  # link to the file already exist


    # Go back to the parent directory and continue.
    os.chdir('..')


