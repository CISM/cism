#!/usr/bin/env python

# This script sets up initial conditions for the MISMIP+ experiment.
# See this paper for details:
#   X. S. Asay-Davis et al. (2016), Experimental design for three interrelated 
#   marine ice sheet and ocean model intercomparison projects: 
#   MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
#   Geosci. Model Devel., 9, 2471-2497, doi: 10.5194/gmd-9-2471-2016.

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


xDomain = 640000.0     # domain x-dimension (m)
yDomain = 80000.0      # domain y-dimension (m)
xCalve  = 640000.      # calving front location (m)
initThickness = 100.   # initial uniform ice thikcness (m)
accum = 0.3            # uniform accumulation rate (m/yr)

restartfreqSpinup = 1000.    # frequency at which restart file is written (yr)


####################################
# Function used later in the code #
####################################


#This function computes the MISMIP+ bed according to Asay-Davis et al. (2016).
def computeBed(x,y):
    x  = x/1.e3     # km
    y  = y/1.e3     # km
    X  = np.size(x)
    Y  = np.size(y)
    b  = np.zeros((X,Y))
    B0 = -150.     # m
    B2 = -728.8    # m
    B4 = 343.91    # m
    B6 = -50.57    # m
    x_bar = 300.   # km
    x_tilde = x/x_bar
    
    dc = 500.      # m
    fc = 4.        # km
    wc = 24.       # km
    Ly = 80.       # km
    
    Bmax = -720.   # m
    B_x  = B0 + B2*x_tilde**2 + B4*x_tilde**4 + B6*x_tilde**6
    B_y  = dc / (1 + np.exp(-2*(y-Ly/2-wc)/fc)) + dc / (1 + np.exp(2*(y-Ly/2+wc)/fc))
    Bsum = B_x + B_y
    B    = np.maximum(Bsum, Bmax)  # B >= Bmax
    return B



########
# Code #
########

# Parse options.
optparser = OptionParser()

optparser.add_option('-c', '--config', dest='configfile',    type='string', default='mismip+.config.template', help="config file template", metavar="FILE")
optparser.add_option('-e', '--exec',   dest='executable',    type='string', default='cism_driver', help="path to the CISM executable")
optparser.add_option('-x', '--expt',   dest='experiment',    type='string', default= 'all',   help="MISMIP+ experiment(s) to set up", metavar="EXPT")
optparser.add_option('-t', '--tstep',  dest='timestep',      type='float',  default= 0.5,     help="time step (yr)",         metavar="TSTEP")
optparser.add_option('-r', '--res',    dest='resolution',    type='int',    default= 2000,    help="grid resolution (m)",    metavar="RES")
optparser.add_option('-v', '--vlevel', dest='vertlevels',    type='int',    default= 3,       help="no. of vertical levels", metavar="VLEVEL")
optparser.add_option('-a', '--approx', dest='approximation', type='string', default= 'DIVA',  help="Stokes approximation (SSA, DIVA, BP)")
optparser.add_option('-b', '--basal',  dest='basalFriction', type='string', default='Schoof', help="basal friction law (Schoof, Tsai, powerlaw)")
optparser.add_option('-y', '--year',   dest='yearsSpinup',   type='int',    default= 20000,   help="length of spinup run (yr)")

optparser.add_option 

for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()

if options.experiment == 'all':
    experiments = ['Spinup', 'Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']
    print 'Setting up all the MISMIP+ experiments'
elif options.experiment == 'allIce':
    experiments = ['Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice1rax', 'Ice1rrx', 'Ice2rax', 'Ice2rrx']
    print 'Run all the MISMIP+ experiments, excluding Spinup'
elif options.experiment == 'Ice1':
    experiments = ['Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx']
    print 'Run the MISMIP+ Ice1 experiments'
elif options.experiment == 'Ice2':
    experiments = ['Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']
    print 'Run the MISMIP+ Ice2 experiments'
elif options.experiment in ['Spinup', 'Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']:
    experiments = [options.experiment]
    print 'Setting up experiment', options.experiment
else:
    sys.exit('Please specify experiment(s) from this list: all, Spinup, Ice0, Ice1r, Ice1ra, Ice1rr, Ice1rax, Ice1rrx, Ice2r, Ice2ra, Ice2rr, Ice2rax, Ice2rrx')

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

print 'MISMIP+ grid resolution (m) =', options.resolution
print 'Number of vertical levels =', nz

# Set number of grid cells in each direction.
# Include a few extra cells in the x direction to handle boundary conditions.
nx = int(xDomain/dx) + 4
ny = int(yDomain/dy)

# Copy the config template to a new master config file.
masterConfigFile = 'mismip+.config'

try:
    shutil.copy(options.configfile, masterConfigFile)
except:
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
if options.basalFriction == 'Schoof':
    which_ho_babc = 11
    print 'Using Schoof basal friction law'
elif options.basalFriction == 'Tsai':
    which_ho_babc = 12
    print 'Using Tsai basal friction law'
elif options.basalFriction == 'powerlaw':
    which_ho_babc = 9
    print 'Using basal friction power law'
else:
    which_ho_babc = 11   # Schoof is default
    print 'Defaulting to Schoof basal friction law'

config.set('ho_options', 'which_ho_babc', which_ho_babc)

yearsSpinup = float(options.yearsSpinup)
config.set('time', 'tend', yearsSpinup)

# Write to the master config file.
with open(masterConfigFile, 'w') as configfile:
    config.write(configfile)

print 'years of spinup =', yearsSpinup
print 'spinup restart frequency =', restartfreqSpinup

# Create the netCDF input file.
try:
    parser = ConfigParser()
    parser.read(options.configfile)
    initfile = parser.get('CF input', 'name')
except:
    sys.exit('Error parsing ' + options.configfile)

print 'Creating input file', initfile
ncfile = Dataset(initfile, 'w')

# Create dimensions.
# Note: (x0,y0) = staggered (velocity) grid.
#       (x1,y1) = unstaggered (scalar) grid.
ncfile.createDimension('time', 1)
ncfile.createDimension('x1',  nx)
ncfile.createDimension('y1',  ny)
ncfile.createDimension('x0',  nx-1)
ncfile.createDimension('y0',  ny-1)
ncfile.createDimension('level',         nz)
ncfile.createDimension('staglevel',     nz-1)
ncfile.createDimension('stagwbndlevel', nz+1)

# Create time and grid variables.
# Note: (x1,y1) are loadable and need to be in the input file.  
#       (x0,y0) are not loadable, but are derived in CISM from (x1,y1). May not be needed.

ncfile.createVariable('time','f4',('time',))[:] = [0]
x1 = ncfile.createVariable('x1','f4',('x1',)) 
y1 = ncfile.createVariable('y1','f4',('y1',))
x0 = ncfile.createVariable('x0','f4',('x0',))
y0 = ncfile.createVariable('y0','f4',('y0',))

# Create 2D input fields
thk  = ncfile.createVariable('thk',  'f4', ('time','y1','x1'))
topg = ncfile.createVariable('topg', 'f4', ('time','y1','x1'))
acab = ncfile.createVariable('acab', 'f4', ('time','y1','x1'))
uvel = ncfile.createVariable('uvel', 'f4', ('time','level','y0','x0'))
vvel = ncfile.createVariable('vvel', 'f4', ('time','level','y0','x0'))
kinbcmask = ncfile.createVariable('kinbcmask', 'i4', ('time','y0','x0'))

# Compute x and y on each grid.
# Note: (1) The x origin is placed at the center of the second cell from the left.
#           This assumes that kinbcmask = 1 at the first vertex from the left.
#           Thus the left edge of the grid has x = -3*dx/2.
#       (2) The y origin is placed at the bottom edge of the CISM grid.
#           The line of central symmetry runs along cell edges at y = 40 km.

x = dx * np.arange(nx,dtype='float32')   # x = 0, dx, 2*dx, etc.
y = dy * np.arange(ny,dtype='float32')   # y = 0, dy, 2*dy, etc.

x1[:] = x[:] - dx         # x1 = -dx, 0, dx, ..., (nx-2)*dx - dx/2  
y1[:] = y[:] + dy/2       # y1 = dy/2, 3*dy/2, ..., (ny-1)*dy - dy/2

x0[:] = x[:-1] - dx/2.    # x0 = -dx/2, dx/2, 3*dx/2, ..., (nx-2)*dx
y0[:] = y[:-1] + dy       # y0 = dy, 2*dy, ..., (ny-1)*dy

# Initialize thickness.
thk[:,:,:] = 0.
for i in range(nx):
    if x1[i] <= xCalve:
        thk[:,:,i] = initThickness


# Set bed topography.
for i in range(nx):
    for j in range(ny):
        topg[:,j,i] = computeBed(x1[i], y1[j])

# Set the surface mass balance.
# Uniform accumulation, but prescribe a large negative rate beyond the calving front.
# WHL - The large negative rate may not be needed, but setting it just in case.
acab[:,:,:] = accum
for i in range(nx):
    if x1[i] > xCalve:
        acab[:,:,i] = -100.   # m/yr

# Set initial velocity to zero.
#WHL- Probably not necessary.
uvel[:,:,:,:] = 0.
vvel[:,:,:,:] = 0.

# Set kinematic velocity mask.
# Where kinbcmask = 1, the velocity is fixed at its initial value.
# Note: Although there is no ice on the RHS of the domain, we need kinbcmask =1 there 
#       to preserve symmetry with the LHS (since east-west BCs are formally periodic).
kinbcmask[:,:,:]  = 0   # initialize to 0 everywhere
kinbcmask[:,:,0]  = 1   # mask out left-most column
kinbcmask[:,:,-1] = 1   # mask out right-most column

ncfile.close()

print 'Experiments:', experiments

# Loop through experiments.
for expt in experiments:

    # For each experiment, make a suitable config file and set up a subdirectory.
    print 'Creating config file for experiment', expt

    # Make a copy of the mismip+Init config file.
    # Below, this copy will be tailored for the chosen MISMIP+ experiment,
    #  without changing the settings used for spin-up.

    newConfigFile = 'mismip+' + expt + '.config'
    print 'Config file for this experiment:', newConfigFile
    shutil.copy(masterConfigFile, newConfigFile)

    # Read the new config file.
    config = ConfigParser()
    config.read(newConfigFile)

    # Experiment-specific settings.
    # Note: The standard experiments are Ice0, Ice1r, Ice1ra, Ice1rr, Ice2r, Ice2ra and Ice2rr.
    #       Experiments Ice1ra, Ice1rr, Ice2ra and Ice2rr are assumed to end at year 200.
    #       Experiments Ice1rax, Ice1rrx, Ice2rax and Ice2rrx are the optional extensions
    #        from year 200 to year 1000.
    # Note: Although these experiments read restart files at startup, they are formally 
    #       cold-start experiments (restart = 0, starting from the 'CF input' file)
    #       rather than restart experiments (restart = 1, starting from the 'CF restart' file).

    if (expt == 'Spinup'):
        tstart      = 0.0
        tend        = yearsSpinup
        inputdir    = '../'
        inputfile   = initfile
        inputslice  = 1
        outputfreq  = min(1000.0, yearsSpinup)
        restartfreq = min(restartfreqSpinup, yearsSpinup)
    elif expt == 'Ice0':
        tstart      = 0.0
        tend        = 100.0
        inputdir    = '../Spinup/'
        inputfile   = 'mismip+Spinup.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice1r':
        config.set('options', 'bmlt_float', 1)
        tstart      = 0.0
        tend        = 100.0
        inputdir    = '../Spinup/'
        inputfile   = 'mismip+Spinup.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice1ra':
        tstart      = 100.0
        tend        = 200.0
        inputdir    = '../Ice1r/'
        inputfile   = 'mismip+Ice1r.restart.nc'
        inputslice  = 1
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice1rax':
        tstart      = 200.0
        tend        = 1000.0
        inputdir    = '../Ice1ra/'
        inputfile   = 'mismip+Ice1ra.restart.nc'
        inputslice  = 1
        outputfreq  = 100.0
        restartfreq = 800.0
    elif expt == 'Ice1rr':
        config.set('options', 'bmlt_float', 1)
        tstart      = 100.0
        tend        = 200.0
        inputdir    = '../Ice1r/'
        inputfile   = 'mismip+Ice1r.restart.nc'
        inputslice  = 1
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice1rrx':
        config.set('options', 'bmlt_float', 1)
        tstart      = 200.0
        tend        = 1000.0
        inputdir    = '../Ice1rr/'
        inputfile   = 'mismip+Ice1rr.restart.nc'
        inputslice  = 1
        outputfreq  = 100.0
        restartfreq = 800.0
    elif expt == 'Ice2r':
        config.set('options', 'bmlt_float', 2)
        tstart      = 0.0
        tend        = 100.0
        inputdir    = '../Spinup/'
        inputfile   = 'mismip+Spinup.restart.nc'
        inputslice  = int(yearsSpinup/restartfreqSpinup)
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice2ra':
        tstart      = 100.0
        tend        = 200.0
        inputdir    = '../Ice2r/'
        inputfile   = 'mismip+Ice2r.restart.nc'
        inputslice  = 1
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice2rax':
        tstart      = 200.0
        tend        = 1000.0
        inputdir    = '../Ice2ra/'
        inputfile   = 'mismip+Ice2ra.restart.nc'
        inputslice  = 1
        outputfreq  = 100.0
        restartfreq = 800.0
    elif expt == 'Ice2rr':
        config.set('options', 'bmlt_float', 2)
        tstart      = 100.0
        tend        = 200.0
        inputdir    = '../Ice2r/'
        inputfile   = 'mismip+Ice2r.restart.nc'
        inputslice  = 1
        outputfreq  = 10.0
        restartfreq = 100.0
    elif expt == 'Ice2rrx':
        config.set('options', 'bmlt_float', 2)
        tstart      = 200.0
        tend        = 1000.0
        inputdir    = '../Ice2rr/'
        inputfile   = 'mismip+Ice2rr.restart.nc'
        inputslice   = 1
        outputfreq  = 100.0
        restartfreq = 800.0


    # Set the time step in the master config file.
    # Set the diagnostic interval to the same value (not necessary, but helpful for debugging).
    # Note: this step is necessary when running at resolution coarser that 4 km as the output files
    #       needs to be written every 10 years to satisfy plotting criteria.
    if expt != 'Spinup':
        config.set('time', 'dt',      min(options.timestep, 2.))
        config.set('time', 'dt_diag', min(options.timestep, 2.))
    else:
        config.set('time', 'dt',      options.timestep)
        config.set('time', 'dt_diag', options.timestep)

    # Set the start and end times.
    config.set('time', 'tstart', tstart)
    config.set('time', 'tend',   tend)

    # Change the default comment.
    comment = 'MISMIP+ experiment ' + expt
    config.set('CF default', 'comment', comment)

    # Set input file and time slice in the section '[CF input]'.
    # Note: This method may not be robust for Spinup runs that start and restart.
    #       For this reason, the script mismip+Run.py makes sure the 'time' entry
    #       in [CF input] corresponds to the final time slice.
    print 'Input file:', inputfile
    config.set('CF input', 'name', inputfile)
    config.set('CF input', 'time', inputslice)

    # Set the output filename in the section '[CF output]'.
    outputfile = 'mismip+' + expt + '.out.nc'
    print 'Output file:', outputfile
    config.set('CF output', 'name',      outputfile)
    config.set('CF output', 'frequency', outputfreq)

    # Set restart info in the section '[CF restart]'.
    # Note: Each experiment (except Spinup) writes only one time slice to a restart file.
    restartfile = 'mismip+' + expt + '.restart.nc'
    print 'Restart file:', restartfile
    config.set('CF restart', 'name',       restartfile)
    config.set('CF restart', 'variables', 'restart')
    config.set('CF restart', 'frequency',  restartfreq)
    config.set('CF restart', 'xtype',     'double')
    config.set('CF restart', 'write_init', False)

    # Write to the new config file.
    with open(newConfigFile, 'w') as configfile:
        config.write(configfile)

    # Create a subdirectory named for the experiment, and stage the run there.
    try:
        os.mkdir(expt)
        print 'Created subdirectory', expt
    except:
        print 'Subdirectory', expt, 'already exists'

    os.chdir(expt)

    # Move the config file from the parent directory to the subdirectory.
    shutil.move('../' + newConfigFile, newConfigFile)
    print 'Created config file', newConfigFile

    # Link to the cism_driver executable in the parent directory.
    try:
        os.symlink('../cism_driver', 'cism_driver')
    except:
        pass   # link to cism_driver already exists

    # Link to the input file in the appropriate directory.
    try:
        os.symlink(inputdir + inputfile, inputfile)
    except:
        pass   # link to the input file already exists

    # Go back to the parent directory and continue.
    os.chdir('..')
