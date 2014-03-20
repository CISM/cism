#!/usr/bin/env python

# File to generate an input netcdf file readable by CISM for the stream test case and run the model.

# Import some modules here
import sys, os, numpy

# =====================================
# parameter values to use in this script and the post-processing plotting script.
# =====================================

analytic_solution = 'raymond'  # can be 'raymond' or 'schoof'
kinflag = 1    # apply kinematic bc (analytic soln) at up/downstream ends
#% kinflag = 0;    %% apply 0 vel bc at up/downstream ends

# Physical parameters
rho = 910.0   # ice density kg/m3
g = -9.81     # gravity m/s2
n = 3         # flow law exponent
A = 1e-16     # flow rate factor in Pa^-3 yr^-1
B = A**(-1.0/n)

# Domain parameters
H = 1000.0       # ice thickness
dsdx = -1.0e-3   # bed (and surface) slope in the x-direction (y-direction is flat)
# ice stream width
w = 50.0e3
W = w / 2
taud = rho * g * H * dsdx

# Schoof solution parameters
m = 1.55  # Schoof exponent
L = 1.4e4  # Schoof yield stress width parameter  ## TODO Steve is this correct?

# -- Functions for analytic solutions --
# Schoof yield stress distribution
def schoof_tau(taud, yy):
  return taud * numpy.absolute( yy / L )**m

# yield stress (for Raymond solution only)
def raymond_tau(x):
  return 5.2e3*numpy.ones(x.shape)


# =====================================
# The actual script to run
# =====================================
if __name__ == '__main__':
    # Parse options
    from optparse import OptionParser
    optparser = OptionParser()
    optparser.add_option("-c", "--config", dest="configfile", type='string', default='stream.config', help="Name of .config file to use to setup and run the test case", metavar="FILE")
    optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel', metavar="NUMPROCS")
    optparser.add_option('-e','--exec',dest='executable',default='./simple_glide',help='Set path to the CISM executable (defaults to "./simple_glide")')
    options, args = optparser.parse_args()


    from netCDF import *
    from math import sqrt
    import ConfigParser


    # =====================================
    # Create a netCDF file according to the information in the config file.
    try:
        configParser = ConfigParser.SafeConfigParser()
        configParser.read(options.configfile)
        nx = int(configParser.get('grid','ewn'))
        ny = int(configParser.get('grid','nsn'))
        nz = int(configParser.get('grid','upn'))
        dx = float(configParser.get('grid','dew'))
        dy = float(configParser.get('grid','dns'))
        filename = configParser.get('CF input', 'name')
    except:
        sys.exit('Error parsing ' + options.configfile)

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

    x1 = dx*numpy.arange(nx,dtype='float32')
    y1 = dx*numpy.arange(ny,dtype='float32')
    x0 = dx/2 + x1[:-1] # staggered grid
    y0 = dy/2 + y1[:-1]

    netCDFfile.createVariable('time','f',('time',))[:] = [0]
    netCDFfile.createVariable('x1','f',('x1',))[:] = x1
    netCDFfile.createVariable('y1','f',('y1',))[:] = y1
    netCDFfile.createVariable('x0','f',('x0',))[:] = x0
    netCDFfile.createVariable('y0','f',('y0',))[:] = y0

    # Calculate values for the required variables.
    thk  = numpy.zeros([1,ny,nx],dtype='float32')
    topg = numpy.zeros([1,ny,nx],dtype='float32')
    tauf = numpy.zeros([1,ny-1,nx-1],dtype='float32')

    # =======================================
    # Calculate input field values
    thk[:] = H  # constant thickness

    for j in range(ny):
      topg[0,j,:] = 1000.0 + dsdx * x1[:]   # sloped bed.  add 1000.0 to stay well above sea level

    if analytic_solution == 'raymond':
        tau0profile = raymond_tau(y0)
    elif analytic_solution == 'schoof':
        tau0profile = schoof_tau(y0)
    else:
        sys.exit("Error: Invalid value for 'analytic_solution'.")
    for i in range(nx-1):
      tauf[0,:,i] = tau0profile
    # put no slip along lateral boundaries - 3 cells thick along north and south
    tauf[0,0:3,:] = 0.7e5
    tauf[0,-3:,:] = 0.7e5

    # =======================================
    # Save the required variables to the netCDF file.
    netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
    netCDFfile.createVariable('tauf','f',('time','y0','x0'))[:] = tauf

    netCDFfile.close()

    # Note: To implement the kinematic bc mask, we will also need to add & populate these variables: 
    #kinbcmask = nc.createVariable('kinbcmask','f',('time','y0','x0'))
    #uvel   = nc.createVariable('uvel',  'f',('time','level','y0','x0'))
    #vvel   = nc.createVariable('vvel',  'f',('time','level','y0','x0'))

    # =======================================
    # Update offset in config file
    offset = -dsdx * dx * nx
    configParser.set('parameters', 'periodic_offset_ew', str(offset))
    configFile = open(options.configfile,'w')
    configParser.write(configFile)
    configFile.close()



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
             mpiexec = 'openmpirun -np '
          elif os.system('which mpirun > /dev/null') == 0:
             mpiexec = 'mpirun -np '
          elif os.system('which aprun > /dev/null') == 0:
             mpiexec = 'aprun -n '
          else:
             sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./simple_glide stream.config')
          runstring = mpiexec + str(options.parallel) + ' ' + options.executable + ' ' + options.configfile
          print 'Executing parallel run with:  ' + runstring + '\n\n'
          os.system(runstring)  # Here is where the parallel run is actually executed!




