#!/usr/bin/env python

# This script reads a CISM output file from a MISMIP experiment, and creates a netCDF output file.
# The input file should include the following:
#      Coordinates: x0, y0, x1, y1, time
#      Global scalars: ivol, ice_mass_above_flotation, iareag
#      2D fields: f_ground, stagthck, usfc, vsfc, ubas, vbas, uvel_mean, vvel_mean
# The f_ground field is not part of the new output file, but is needed to identify
# vertices lying on the grounding line (GL).


import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from optparse import OptionParser


###############################
# Constants used in this code #
###############################

#Set ice density prescribed for MISMIP3d.
rhoi  = 900.      # kg/m^3
model = '_cism'   # file naming extension


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



########
# Code #
########


# Parse options.
parser = OptionParser()

parser.add_option('-x', '--expt',dest='experiment',type='string', default='all',   help='MISMIP experiment(s) to run', metavar="EXPT")
parser.add_option('-f', '--file',dest='filename',  type='string',                  help='CISM output file from MISMIP3d run', metavar="FILE")
parser.add_option('-s', '--stat',dest='StatChoice',type='string',default='advance',help='MISMIP experiment set to run', metavar="EXPT")
parser.add_option('--bed', dest='bedtopo', type='string', default ='linear',help='bed topography, linear or poly', metavar='BEDTOPO')


for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"

options, args = parser.parse_args()

if options.bedtopo == 'linear':
    As             = AsLinear
    AsAdvance      = AsAdvanceLinear
    AsRetreat      = AsRetreatLinear
    Astatus        = AstatusLinear
    AstatusAdvance = AstatusAdvanceLinear
    AstatusRetreat = AstatusRetreatLinear
elif options.bedtopo == 'poly':
    As             = AsPoly
    AsAdvance      = AsAdvancePoly
    AsRetreat      = AsRetreatPoly
    Astatus        = AstatusPoly
    AstatusAdvance = AstatusAdvancePoly
    AstatusRetreat = AstatusRetreatPoly
else:
    sys.exit('Please specify bed type from this list: linear, poly')


if options.experiment == 'all':
    experiments = As
    Astat       = Astatus
    print 'Writing all the MISMIP experiments'
elif options.experiment == 'advance':
    experiments = AsAdvance
    Astat       = AstatusAdvance
    print 'Writing advance experiments'
elif options.experiment == 'retreat':
    experiments = AsRetreat
    Astat       = AstatusRetreat
    print 'Writing retreat experiments'
else:
    sys.exit('Please specify experiment(s) from this list: all, advance, retreat')



# Loop through A values.
bedType   = options.bedtopo
countStat = -1               # counter to access status matrix. Accounting for zero-arrays indexing.

for expt in experiments:
    
    countStat = countStat + 1
    stat      = Astat[countStat]
    
    # Change to bed type directory.
    os.chdir(bedType)

    # Change to Advance or Retreat directory.
    os.chdir(stat)

    # Change to the subdirectory for this experiment.
    os.chdir(expt)

    if options.filename:
        file = options.filename
    else:
        file = 'mismip_' + expt + '.out.nc'

    print 'Creating a MISMIP grounding-line file for experiment', expt
    print 'Attempting to read CISM file', file

    # Open the CISM output file, get needed dimensions.
    # Note: (x0,y0) are dimensions of the staggered (velocity) grid
    #       (x1,y1) are dimensions of the unstaggered (scalar) grid

    try:
        cismfile = Dataset(file,'r')
    except:
        sys.exit('Error: Unable to open CISM file')

    try:
        nTime = len(cismfile.dimensions['time'])
        nx    = len(cismfile.dimensions['x0'])
        ny    = len(cismfile.dimensions['y0'])
    except:
        sys.exit('Error: The CISM file is missing needed dimensions')

   # Initialize some variables and arrays.
   # Read in some fields needed to compute diagnostics similar to the one of MISMIP+.

    print 'Reading in CISM variables...'

    try:
        # These array names are somewhat arbitrary.  Sometime I have used CISM names;
        # sometimes I have used a different name that is closer (but not identical)
        # to the MISMIP3d field names.

        # Read horizontal and time coordinates.
        # Note: (x0,yo) are vertex coordinates; (x1,y1) are cell center coordinates.
        x0 = cismfile.variables['x0'][:]
        y0 = cismfile.variables['y0'][:]
        x1 = cismfile.variables['x1'][:]
        y1 = cismfile.variables['y1'][:]
        t  = cismfile.variables['time'][:]

        # Read scalars (function of time only).
        ivol                  = cismfile.variables['ivol'][:]
        imass_above_flotation = cismfile.variables['imass_above_flotation'][:]
        iareag                = cismfile.variables['iareag'][:]

        # Read 2D fields (functions of x, y and time).
        # Note: All these fields are located at vertices.
        f_ground     = cismfile.variables['f_ground'][:,:,:]
        iceThickness = cismfile.variables['thk'][:,:,:]
        uSurface     = cismfile.variables['usfc'][:,:,:]
        vSurface     = cismfile.variables['vsfc'][:,:,:]
        uBase        = cismfile.variables['ubas'][:,:,:]
        vBase        = cismfile.variables['vbas'][:,:,:]
        uMean        = cismfile.variables['uvel_mean'][:,:,:]
        vMean        = cismfile.variables['vvel_mean'][:,:,:]

    except:
        sys.exit('Error: The output file is missing needed fields.')

    # Create the GL output file.

    outfilename = expt + model + '.nc'
    ncfile = Dataset(outfilename, 'w')
    print 'Created output file', outfilename

    # Set dimensions.
    glptdim = ncfile.createDimension('nPointGL', size = None)
    timedim = ncfile.createDimension('nTime', size = nTime)

    # Create variables.
    xGL  = ncfile.createVariable('xGL',  'f4', ('nPointGL', 'nTime'))
    yGL  = ncfile.createVariable('yGL',  'f4', ('nPointGL', 'nTime'))
    time = ncfile.createVariable('time', 'f4', ('nTime'))

    iceVolume    = ncfile.createVariable('iceVolume',    'f4', ('nTime'))
    iceVAF       = ncfile.createVariable('iceVAF',       'f4', ('nTime'))
    groundedArea = ncfile.createVariable('groundedArea', 'f4', ('nTime'))

    iceThicknessGL = ncfile.createVariable('iceThicknessGL', 'f4', ('nPointGL', 'nTime'))
    uSurfaceGL     = ncfile.createVariable('uSurfaceGL',     'f4', ('nPointGL', 'nTime'))
    vSurfaceGL     = ncfile.createVariable('vSurfaceGL',     'f4', ('nPointGL', 'nTime'))
    uBaseGL        = ncfile.createVariable('uBaseGL',        'f4', ('nPointGL', 'nTime'))
    vBaseGL        = ncfile.createVariable('vBaseGL',        'f4', ('nPointGL', 'nTime'))
    uMeanGL        = ncfile.createVariable('uMeanGL',        'f4', ('nPointGL', 'nTime'))
    vMeanGL        = ncfile.createVariable('vMeanGL',        'f4', ('nPointGL', 'nTime'))

    # Loop over time slices and fill variables.
    print 'Adding grounding-line variables to output file...'

    for iTime in range(nTime):

        print '   Time slice:', iTime

        # Add the scalar data for this time slice.
        time[iTime]         = t[iTime]
        iceVolume[iTime]    = ivol[iTime]
        iceVAF[iTime]       = imass_above_flotation[iTime]/rhoi
        groundedArea[iTime] = iareag[iTime]

        # Loop over the horizontal grid and identify GL points.
        # For each GL point, add the desired data to the arrays.
        eps = 1.0e-6
        nGL = 0

        # A cell hosts a grounding line if at least one neighboring point has a zero f_ground value.
        for i in range(1,nx-1):
            if (f_ground[iTime,3,i] > 0) and (f_ground[iTime,3,i+1] == 0):
                    nGL = nGL + 1
                    m = nGL - 1  # indexing starts at 0
                    xGL[m,iTime] = x0[i]
                    
                    iceThicknessGL[m,iTime] = iceThickness[iTime,3,i]
                    uSurfaceGL[m,iTime]     = uSurface[iTime,3,i]
                    vSurfaceGL[m,iTime]     = vSurface[iTime,3,i]
                    uBaseGL[m,iTime]        = uBase[iTime,3,i]
                    vBaseGL[m,iTime]        = vBase[iTime,3,i]
                    uMeanGL[m,iTime]        = uMean[iTime,3,i]
                    vMeanGL[m,iTime]        = vMean[iTime,3,i]



    ncfile.close()

    # Change to the parent directory.
    os.chdir('../../..')
