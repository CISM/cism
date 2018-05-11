#!/usr/bin/env python

# This script sets up initial conditions for the MISMIP (1d) experiment.
# See this paper for details:
# Pattyn et al., Results of the Marine Ice Sheet Model Intercomparison Project, MISMIP,
#    The Cryosphere, 6, 573-588, 2012, doi:10.5194/tc-6-573-2012.
#
# Note: This experiment is meant to analyse flowline models.
#       There is no perturbation whatsoever in the y direction.
#       Therefore the domain is limited to 6 grid cells in the y direction.

import sys, os
import shutil
import fileinput
import numpy as np
from netCDF4 import Dataset
from ConfigParser import ConfigParser
from optparse import OptionParser


###############################
# Constants used in this code #
###############################

accum = 0.3    # uniform accumulation (m/yr)


#### Linear bed specific ####

# A-values used in the linear bed experiment.
AsLinear = ['14.638e-17','6.7941e-17','3.1536e-17',
            '14.638e-18','6.7941e-18','3.1536e-18',
            '14.638e-19','6.7941e-19','3.1536e-19',
            '6.7941e-19','14.638e-19',
            '3.1536e-18','6.7941e-18','14.638e-18',
            '3.1536e-17','6.7941e-17','14.638e-17']

# A-values used in the linear bed and advance experiment.
AsAdvanceLinear = ['14.638e-17','6.7941e-17','3.1536e-17',
                   '14.638e-18','6.7941e-18','3.1536e-18',
                   '14.638e-19','6.7941e-19','3.1536e-19']

# A-values used in the linear bed and retreat experiment.
AsRetreatLinear = list(reversed(AsAdvanceLinear[0:-1]))


# Status of the linear bed experiment.
AstatusLinear = ['advance','advance','advance',
                 'advance','advance','advance',
                 'advance','advance','advance',
                 'retreat','retreat',
                 'retreat','retreat','retreat',
                 'retreat','retreat','retreat']

# Status of the linear bed and advance experiment.
AstatusAdvanceLinear = ['advance','advance','advance',
                        'advance','advance','advance',
                        'advance','advance','advance']

# Status of the linear bed and retreat experiment.
AstatusRetreatLinear = ['retreat','retreat',
                        'retreat','retreat','retreat',
                        'retreat','retreat','retreat']



#### Poly bed specific ####


# A-values used in the poly bed experiment.
AsPoly = ['9.4608e-18','7.8840e-18','6.3072e-18',
          '4.7304e-18','3.1536e-18','1.5768e-18',
          '7.8840e-19',
          '1.5768e-18','3.1536e-18','4.7304e-18',
          '6.3072e-18','7.8840e-18','9.4608e-18']

# A-values used in the poly bed and advance experiment.
AsAdvancePoly = ['9.4608e-18','7.8840e-18','6.3072e-18',
                 '4.7304e-18','3.1536e-18','1.5768e-18',
                 '7.8840e-19']

# A-values used in the poly bed and retreat experiment.
AsRetreatPoly = list(reversed(AsAdvancePoly[0:-1]))


# Status of the poly bed experiment.
AstatusPoly = ['advance','advance','advance',
               'advance','advance','advance',
               'advance',
               'retreat','retreat','retreat',
               'retreat','retreat','retreat']

# Status of the poly bed and advance experiment.
AstatusAdvancePoly = ['advance','advance','advance',
                      'advance','advance','advance',
                      'advance']

# Status of the poly bed and retreat experiment.
AstatusRetreatPoly = ['retreat','retreat','retreat',
                      'retreat','retreat','retreat']


# Final times prescribed for the poly bed experiment.
AsTimePoly = ['30000','15000','15000',
              '15000','15000','30000',
              '30000',
              '15000','15000','30000',
              '30000','30000','15000']

# Final times prescribed for the poly bed experiment.
AsTimeAdvancePoly = ['30000','15000','15000',
                     '15000','15000','30000',
                     '30000']

# Final times prescribed for the poly bed experiment.
AsTimeRetreatPoly = ['15000','15000','30000',
                     '30000','30000','15000']



####################################
# Function used later in the code #
####################################


# The following function returns the linear bed topography as in Pattyn et al. (2012).
def computeBedLinear(x):
    # Input expected in km.
    
    schoofx = 750.    # scaling factor in km
    slope   = -778.5  # m
    b0      = 720.    # m
    
    eps_b  = 1e-10
    abs_x  = np.sqrt(x**2 + eps_b**2)  # km (smoothing for ice divide BC requirements)
    xprime = abs_x/schoofx             # unitless
    b      = b0 + slope*xprime         # m
    
    return b

# The following function returns the polynomial bed topography as in Pattyn et al. (2012).
def computeBedPoly(x):
    # Input x expected in km.
    
    schoofx = 750.    # scaling factor in km
    x  = x/schoofx    # unitless
    
    b0 = 729.         # m
    b2 = -2184.8      # m
    b4 = 1031.72      # m
    b6 = -151.72      # m
    
    b  = b0 + b2*x**2 + b4*x**4 + b6*x**6    # m
    
    return b

# The following function returns the bed topography as in Pattyn et al. (2012) based
# on a choice between 'linear' and 'polynomial'
def computeBed(x, bedType):
    # Input x expected in km.
    
    if bedType == 'linear':
        b = computeBedLinear(x)
    elif bedType == 'poly':
        b = computeBedPoly(x)
    else:
        sys.exit('Please specify bed topography from these options: linear or poly.')
                 
    return b


# The following function computes a uniform initial ice thickness for the experiment.
def computeHGuessSlab():
    
    H = 100.    # m
    return H


# The following function computes an initial ice thickness based on semi-analyitc
# solution of Schoof2007.
def computeHGuessSemiAnalytic(xH, xu, bedType):
     # Input xH and xu expected in m.

     sPerY = 365.0*24.*3600. # number of seconds per year
     n     = 3.
     rhoi  = 900.            # kg/m^3
     rhow  = 1000.           # kg/m^3
     g     = 9.8             # m/s^2
     delta = 1 - rhoi/rhow   # unitless
     a     = 0.3/sPerY       # converting accumulation to m/s
     C     = 7.624e6         # Pa(s/m)^(1/3)
     A     = 4.6416e-24      # Pa^-3/s
     xg    = 900000.0        # m

     deltaX = np.abs(xH[1]-xH[0])  # m

     # Determining the index of GL relative to domain.
     glIndex = np.argmin(np.abs(xH[:]-xg))
     xg = xH[glIndex]
    
     # Computing bed depth and ice thickness at the grounding line.
     # Note: the sign convention here is different than in Schoof2007.
     bxg = -computeBed(xg*1.e-3, bedType)
     Hxg = bxg/(1.0-delta)
     
     if(Hxg <= 0.):
         raise ValueError("H(xg) <= 0. Cannot produce HGuess")
    
     b   = -computeBed(xH[:]*1e-3, bedType)
     uxg = a*xg/Hxg
     xH  = xH[:]

     # Computing the analytic solution in ice shelf.
     c1 = rhoi*delta*g*a*A**(1./n)/4
     c1 = c1**n
     operand = np.maximum(c1*(xH[:]**(n+1) - xg**(n+1)),0.0) + uxg**(n+1)
     uShelf  = (operand)**(1/(n+1))
     HShelf  = a*xH[:]/uShelf[:]

     H = HShelf.copy()
     u = uShelf.copy()
    
     # Assume balance between taub and taud in the sheet (SIA) and steady state (u=a*x/H),
     # and iteratively solve the ODE for H
     Hip1 = H[glIndex]

     deltaB = b[5]-b[4]
     for xIndex in range(glIndex-1,-1,-1):
         deltaB    = (b[xIndex+1]-b[xIndex])/deltaX
         Hi        = Hip1+deltaX*(deltaB + C/(rhoi*g)*np.abs(a*xH[xIndex])**(1./n)/(Hip1**(1./n+1)))
         Hip1      = Hi
         H[xIndex] = Hi

     H[0] = H[1]   # Enforcing ice divide boundary condition.

     return H


# This function returns the initial thickness profile for the experiment. The choices are between
# an initial uniform slab and semi-analytical solution from Schoof2007 given an initial ad-hoc GL position.
def computeHGuess(xH, xu, bedType, HInitType):
    # Input xH and xu expected in m.
    
    if HinitType == 'slab':
        H = computeHGuessSlab()
    elif HinitType == 'analytic':
        H = computeHGuessSemiAnalytic(xH, xu, bedType)
    else:
        sys.exit('Please specify initial profile from these options: slab or analytic.')

    return H



#############
# Main code #
#############


# Parse options.
optparser = OptionParser()

optparser.add_option('-c', '--config',dest='configfile',   type='string',default='mismip.config.template', help='config file name for setting up the MISMIP experiment', metavar='FILE')
optparser.add_option('-e', '--exec',  dest='executable',default='cism_driver',help='Set path to the CISM executable', metavar='EXECUTABLE')
optparser.add_option('-x', '--expt',  dest='experiment',   type='string',default = 'all', help='MISMIP experiment(s) to set up', metavar='EXPT')
optparser.add_option('-t', '--tstep', dest='timestep',     type='float', default = 1,     help='time step (yr)', metavar='TSTEP')
optparser.add_option('-r', '--res',   dest='resolution',   type='int',   default = 2000,  help='grid resolution (m)', metavar='RES')
optparser.add_option('-v', '--vlevel',dest='vertlevels',   type='int',   default = 3,     help='no. of vertical levels', metavar='VLEVEL')
optparser.add_option('-a', '--approx',dest='approximation',type='string',default = 'DIVA',help='Stokes approximation (SSA, DIVA, BP)', metavar='APPROXIMATION')
optparser.add_option('-b', '--basal', dest='basalFriction', type='string', default='powerlaw', help='Basal friction law (powerlaw, schoof)', metavar='BASALFRICTION')
optparser.add_option('-y', '--year',  dest='yearsSpinup', type='int', default = 30000, help='Length of Spinup run (yr)', metavar='YEARSPINUP')
optparser.add_option('--bed',  dest='bedtopo',  type='string',default ='linear',help='bed topography, linear or poly', metavar='BEDTOPO')
optparser.add_option('--yrun', dest='yearsRun', type='int',   default ='20000', help='run length between 2 experiments', metavar='YEARSRUN')
optparser.add_option('--hinit',dest='initThick',type='string',default ='slab',  help='experiment initial thickness profile', metavar='INITTHICK')

optparser.add_option

for option in optparser.option_list:
    if option.default != ('NO', 'DEFAULT'):
        option.help += (' ' if option.help else '') + '[default: %default]'
options, args = optparser.parse_args()


if options.bedtopo == 'linear':
    As             = AsLinear
    AsAdvance      = AsAdvanceLinear
    AsRetreat      = AsRetreatLinear
    Astatus        = AstatusLinear
    AstatusAdvance = AstatusAdvanceLinear
    AstatusRetreat = AstatusRetreatLinear
    
    xDomain      = 1840000.0    # domain x-dimension (m)
    marine_limit = -1140        # configuration to chop off ice passed this depth
elif options.bedtopo == 'poly':
    As             = AsPoly
    AsAdvance      = AsAdvancePoly
    AsRetreat      = AsRetreatPoly
    Astatus        = AstatusPoly
    AstatusAdvance = AstatusAdvancePoly
    AstatusRetreat = AstatusRetreatPoly

    xDomain      = 1520000.0    # domain x-dimension (m)
    marine_limit = -1100        # configuration to chop off ice passed this depth
else:
    sys.exit('Please specify bed type from this list: linear, poly.')


if options.experiment == 'all':
    experiments = As
    AsTime      = AsTimePoly
    Astat       = Astatus
    print 'Setting up all the experiments'
elif options.experiment == 'advance':
    experiments = AsAdvance
    AsTime      = AsTimeAdvancePoly
    Astat       = AstatusAdvance
    print 'Setting up advance experiments'
elif options.experiment == 'retreat':
    experiments = AsRetreat
    AsTime      = AsTimeRetreatPoly
    Astat       = AstatusRetreat
    print 'Setting up retreat experiments'
else:
    sys.exit('Please specify experiment(s) from this list: all, advance, retreat.')


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

print 'MISMIP grid resolution (m) =', options.resolution
print 'Number of vertical levels =', nz

# Note: This is a streamline experiment with no y-direction variation in bed or forces.
#       For this reason we can limit the domain in y-direction by a fixed amount of grid cells.
yDomain = dy*5

# Set number of grid cells in each direction.
# Include a few extra cells in the x direction to handle boundary conditions.
nx = int(xDomain/dx) + 4
ny = int(yDomain/dy)

# Copy the config template to a new master config file.
masterConfigFile = 'mismip.config'

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

config.set('time', 'dt',      options.timestep)
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
    which_ho_babc = 9   # powerlaw is default
    print 'Defaulting to powerlaw basal friction law'

config.set('ho_options', 'which_ho_babc', which_ho_babc)

# Config setting related to spin up time.
yearsSpinup = float(options.yearsSpinup)
config.set('time', 'tend', yearsSpinup)

# Config setting related to presence of ice.
config.set('parameters','marine_limit',marine_limit)

# Write to the master config file.
with open(masterConfigFile, 'w') as configfile:
    config.write(configfile)


print 'years of Spinup experiment =', yearsSpinup
restartfreqSpinup = min(1000.0, options.yearsSpinup)    # can be changed by the user if needed
print 'Spinup restart frequency =', restartfreqSpinup


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
x1 = ncfile.createVariable('x1','f4',('x1',))
y1 = ncfile.createVariable('y1','f4',('y1',))
x0 = ncfile.createVariable('x0','f4',('x0',))
y0 = ncfile.createVariable('y0','f4',('y0',))

# Create 2D input fields.
thk  = ncfile.createVariable('thk',  'f4', ('time','y1','x1'))
topg = ncfile.createVariable('topg', 'f4', ('time','y1','x1'))
acab = ncfile.createVariable('acab', 'f4', ('time','y1','x1'))
uvel = ncfile.createVariable('uvel', 'f4', ('time','level','y0','x0'))
vvel = ncfile.createVariable('vvel', 'f4', ('time','level','y0','x0'))
kinbcmask = ncfile.createVariable('kinbcmask', 'i4', ('time','y0','x0'))  # kinematic BC mask

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
bedType = options.bedtopo
print 'Computing ' + bedType + ' bed'
for i in range(nx):
    topg[:,:,i] = computeBed(x1[i]/1.e3, bedType) # x1 is in [m] and we need input in [km]


# Creating 'linear' or 'poly' bedtype directory.
try:
    os.mkdir(bedType)
    print 'Created subdirectory', bedType
except OSError:
    print 'Subdirectory', bedType, 'already exists'


# Set initial ice thickness.
HinitType = options.initThick
initThickness = computeHGuess(x1, x0, bedType, HinitType)
for j in range(ny):
    thk[0,j,:] = initThickness

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


ncfile.close()


# Loop through A values.
AprevString = ''    # string used in linking input file of the experiments.
countTime   = -1    # counter to access end time matrix for poly bed. Accounting for zero-array indexing
countStat   = -1    # counter to access status matrix. Accounting zero-array indexing

for expt in experiments:
    
    countTime = countTime + 1
    countStat = countStat+1
    stat      = Astat[countStat]

    # Change to bed type directory.
    os.chdir(bedType)
    
    # Create advance and/or retreat directory if needed.
    try:
        os.mkdir(stat)
        print 'Created subdirectory', expt
    except OSError:
        print 'Subdirectory', expt, 'already exists'

    # Change to advance or retreat directory.
    os.chdir(stat)

    # For each experiment, make a suitable config file and set up a subdirectory.
    print 'Creating config file for experiment A value', expt
    
    # Make a copy of the mismip config file.
    # Below, this copy will be tailored for the chosen MISMIP experiment,
    # without changing the settings used for the Spinup experiment.
    
    newConfigFile = 'mismip_' + expt + '.config'
    print 'Config file for this experiment:', newConfigFile
    shutil.copy('../../' + masterConfigFile, newConfigFile)
    
    # Read the new config file.
    config = ConfigParser()
    config.read(newConfigFile)

    # Experiment-specific settings.
    if bedType == 'linear':
        if (expt == As[0]) and (stat == 'advance'):
            tstart = 0.0
            tend   = yearsSpinup
            inputdir    = '../../../'
            inputfile   = initfile
            inputslice  = 1
            outputfreq  = min(1000.0, restartfreqSpinup)
            restartfreq = restartfreqSpinup
        elif (expt == AsRetreat[0]) and (stat == 'retreat'):
            tstart = 0.0
            tend   = float(options.yearsRun)
            inputdir    = '../../advance/' + AprevString + '/'
            inputfile   = 'mismip_' + AprevString + '.restart.nc'
            inputslice  = 1
            outputfreq  = min(1000.0, restartfreqSpinup)
            restartfreq = restartfreqSpinup
        else:
            tstart = 0.0
            tend   = float(options.yearsRun)
            inputdir    = '../' + AprevString + '/'
            inputfile   = 'mismip_' + AprevString + '.restart.nc'
            inputslice  = 1
            outputfreq  = 1000.0
            restartfreq = 1000.0
    elif bedType == 'poly':
        if (expt == As[0]) and (stat == 'advance'):
            tstart = 0.0
            tend   = AsTime[countTime]
            inputdir    = '../../../'
            inputfile   = initfile
            inputslice  = 1
            outputfreq  = min(1000.0, restartfreqSpinup)
            restartfreq = restartfreqSpinup
        elif (expt == AsRetreat[0]) and (stat == 'retreat'):
            tstart = 0.0
            tend   = AsTime[countTime]
            inputdir    = '../../advance' + AprevString + '/'
            inputfile   = 'mismip_' + AprevString + '.restart.nc'
            inputslice  = 1
            outputfreq  = min(1000.0, restartfreqSpinup)
            restartfreq = restartfreqSpinup
        else:
            tstart = 0.0
            tend   = AsTime[countTime]
            inputdir    = '../' + AprevString + '/'
            inputfile   = 'mismip_' + AprevString + '.restart.nc'
            inputslice  = 1
            outputfreq  = 1000.0
            restartfreq = 1000.0
    else:
        print('This should not be an option by now.')
        sys.exit('Exiting the run.')

    # Set the start and end times.
    config.set('time', 'tstart', tstart)
    config.set('time', 'tend',   tend)

    # Set the default flwa value.
    config.set('parameters', 'default_flwa', float(expt))

    # Change the default comment.
    comment = 'MISMIP experiment ' + expt
    config.set('CF default', 'comment', comment)
    
    # Set input file and time slice.
    config.set('CF input', 'name', inputfile)
    config.set('CF input', 'time', inputslice)

    # Set the output filename in the section [CF output].
    outfilename = 'mismip_' + expt + '.out.nc'
    print 'Output file:', outfilename
    config.set('CF output', 'name',      outfilename)
    config.set('CF output', 'frequency', outputfreq)

    # Set restart info.  This should be in the section called '[CF output]'.
    # Note: Each experiment (except Stnd) writes out only one time slice to a restart file.
    restartfilename = 'mismip_' + expt + '.restart.nc'
    print 'Restart file:', restartfilename
    config.set('CF restart', 'name',       restartfilename)
    config.set('CF restart', 'variables',  'restart')
    config.set('CF restart', 'frequency',  restartfreq)
    config.set('CF restart', 'start',      tstart + restartfreq)
    config.set('CF restart', 'xtype',      'double')

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
        os.symlink('../../../cism_driver', 'cism_driver')
    except OSError:
        pass   # link to cism_driver already exists

    # Link to the input file in the appropriate directory.
    try:
        os.symlink(inputdir + inputfile, inputfile)
    except OSError:
        pass   # link to the input file already exists


    # Updating the previous values of Aprev for next experiment setup.
    AprevString = expt


    # Go back to the parent directory and continue.
    os.chdir('../../..')
