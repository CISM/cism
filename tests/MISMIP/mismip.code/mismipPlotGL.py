#!/usr/bin/env python

# This script plots the grounding line position for the MISMIP experiments.
# This script requires the user to have run the python script "mismipWriteGL.py".


from netCDF4 import Dataset
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import sys, os



###############################
# Constants used in this code #
###############################

model = '_cism'          # file naming extension
sPerY = 365.0*24.*3600.  # number of second in a year

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
    bx     = slope/schoofx*x/abs_x     # m/km
    
    return (b,bx)

# The following function returns the polynomial bed topography as in Pattyn et al. (2012).
def computeBedPoly(x):
    # Input expected in km.
    
    schoofx = 750.    # scaling factor in km
    x  = x/schoofx    # unitless
    
    b0 = 729.         # m
    b2 = -2184.8      # m
    b4 = 1031.72      # m
    b6 = -151.72      # m
    
    b  = b0 + b2*x**2 + b4*x**4 + b6*x**6            # m
    bx = (2*b2*x + 4*b4*x**3 + 6*b6*x**5)/schoofx    # m/km
    
    return (b,bx)


# The following function returns the bed topography as in Pattyn et al. (2012) based
# on the choice between 'linear' and 'polynomial'
def computeBed(x, bedType):
    # Input expected in km.
    
    if bedType == 'linear':
        (b,bx) = computeBedLinear(x)
    elif bedType == 'poly':
        (b,bx) = computeBedPoly(x)
    else:
        sys.exit('Please specify bed topography from these options: linear or poly.')
    
    return (b,bx)


# The following function returns the semi-analytic solution from Schoof2007 model A.
# Note: the units used in Schoof2007 are not the same as the one in CISM.
def xgSemianalytic(bedType):
    
    sPerY = 365.0*24.*3600. # number of seconds per year
    n     = 3.
    rhoi  = 900.            # kg/m^3
    rhow  = 1000.           # kg/m^3
    g     = 9.8             # m/s^2
    delta = 1 - rhoi/rhow   # unitless
    a     = 0.3/sPerY       # converting accumulation to m/s
    C     = 7.624e6         # Pa(s/m)^(1/3)

    if bedType == 'linear':
        xx = np.linspace(1.e6,1.8e6,400)
    else:
        xx = np.linspace(0.7e6,1.5e6,400)


    AmodelA = np.zeros(len(xx))
    for k in range(0,len(xx)):
        x = xx[k]

        (b,bx) = computeBed(x*1.e-3, bedType)
        bx     = bx*1.e-3    # needs to be in m/m for unit compliance with schoof2007

        # Schoof model A as A = f(x).
        h = -b/(1-delta)
        q = x*a
        u = q/h

        num =  -a + u*(-bx-C/(rhoi*g*h)*np.abs(q)**(1./n-1)*q/h**(1./n))
        den = (rhoi*g*delta/4)**n*h**(n+1)
        AmodelA[k] = -num/den

    return (AmodelA,xx)




########
# Code #
########


# Parse options.
optparser = OptionParser()

optparser.add_option('-x', '--expt',     dest='experiment', type='string', default = 'all', help='MISMIP experiment set to run', metavar="EXPT")
optparser.add_option('-s', '--stat',    dest='StatChoice',type='string',default='advance', help='MISMIP experiment set to run', metavar="EXPT")
optparser.add_option('--bed', dest='bedtopo', type='string', default ='linear',help='bed topography, linear or poly', metavar='BEDTOPO')

for option in optparser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
options, args = optparser.parse_args()


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
    print 'Plotting all the MISMIP experiments'
elif options.experiment == 'advance':
    experiments = AsAdvance
    Astat       = AstatusAdvance
    print 'Plotting advance experiments'
elif options.experiment == 'retreat':
    experiments = AsRetreat
    Astat       = AstatusRetreat
    print 'Plotting retreat experiments'
else:
    sys.exit('Please specify experiment(s) from this list: all, advance, retreat')



# Loop through A values.
count   = -1
Aval    = np.zeros(len(experiments))    # initialize A-values storing array
xgval   = np.zeros(len(experiments))    # initialize GL storing array
bedType = options.bedtopo               # type of bed topography

for expt in experiments:
    count = count + 1
    
    Aval[count]  = float(expt)
    stat = Astat[count]

    # Change to bed type directory.
    os.chdir(bedType)

    # Change to Advance or Retreat directory.
    os.chdir(stat)

    # Change to the subdirectory for this experiment.
    os.chdir(expt)

    # Read the file and extract grounding line position.
    try:
       file  = expt + model + '.nc'
       ncid  = Dataset(file, 'r')
       xgval[count] = ncid.variables["xGL"][-1][-1]
       ncid.close()
    except:
       print 'Results for experiment',stat,'and',expt,'is not available'
       

    # Switch back to original directory.
    os.chdir('../../../')


# Obtain the semi-analytic solution from Schoof 2007.
# Note: we compute invAnal (1/A) because it will be used to display the results.
#       It has a better esthetic this way.
(Aanal, xanal) = xgSemianalytic(bedType)
invAanal = 1./Aanal     # computing 1/A of semi-analytic solution

# Compute 1/A of values used for simulation and converting to the
# same units as in PAttyn et al. 2012.
invA = 1./(Aval/sPerY)
invAadv = invA[0:len(AsAdvance)]     # 1/A used in advance experiments
invAret = list(reversed(invAadv))    # 1/A used in retreat experiments
xgadv   = xgval[0:len(AsAdvance)]    # GL of advanced experiments
xgret   = xgval[len(AsAdvance)-1:]   # GL of retreat experiments

# Adjust plot properties depending on the bed type.
if bedType == 'linear':
    ymin   = 1000.0
    ymax   = 1800.0
    nytick = 9
else:
    ymin   = 700.0
    ymax   = 1500.0
    nytick = 9


# Plot the figure displaying the grounding line location.
# Note: the figure will display one subplot for the advance experiments
#       and another one for the retreat experiments.

# Set figure size.
plt.figure(figsize=(7, 7))

# Advance experiment display.
plt.subplot(121)

# Plot semi-analytical solution from Schoof 2007.
plt.semilogx(invAanal, xanal*1.e-3, color='black', label='analytic',linewidth=1)

# Plot simulation results.
plt.semilogx(invAadv,xgadv*1.e-3, '+', ms=10, mfc='red',mec='red',label='simulation advance')

# Turn on logarythmic grid display.
plt.grid(True,which="both",ls="-")

# Set limit display for x and y axis.
plt.xlim((invAadv[0], invAadv[-1]))
plt.ylim((ymin,ymax))

# Add labels for x and y axis and their position and font size.
plt.xlabel("1/A (Pa$^3$s)", position=(1,0), size=12)
plt.ylabel("xg(km)", size=12)

# Set the legend location to upper left.
plt.legend(loc=2)


# Retreat experiment display.
plt.subplot(122)

# Plot semi-analytical solution from Schoof 2007 with reversed x-axis.
plt.semilogx(invAanal, xanal*1.e-3, color='black', label='analytic',linewidth=1)
plt.gca().invert_xaxis()

# Plot simulation results with reversed x-axis.
plt.semilogx(invAret,xgret*1.e-3, '+', ms=10, mfc='blue',mec='blue',label='simulation retreat')
plt.gca().invert_xaxis()

# Turn on logarithmic grid display.
plt.grid(True,which="both",ls="-")

# Set limit display for x and y axis.
plt.xlim((invAret[0], invAret[-1]))
plt.ylim((ymin,ymax))

# Hack to display the gridding of y-axis using display from subplot 1.
plt.yticks(np.linspace(ymin,ymax,nytick)," ")

# Set the legend location to upper right.
plt.legend(loc=1)

# Suppress the vertical white spacing between the 2 subplots.
# This is how they look attached to one another.
plt.subplots_adjust(wspace=0)

# Add a title to the figure based on its location on subplot 2.
plt.title('Grounding line position', position=(0,1))

# Save the figure.
if bedType == 'linear':
    plt.savefig("mismipPlotGLLinearBed.pdf")
elif bedType == 'poly':
    plt.savefig("mismipPlotGLPolyBed.pdf")
else:
    print('Saving the figure with a randome name')
    plt.savefig("mismipPlotGLRandomBed.pdf")
