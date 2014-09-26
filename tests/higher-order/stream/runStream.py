#!/usr/bin/env python

# File to generate an input netcdf file readable by CISM for the stream test case and run the model.

# Import some modules here
import numpy

# =====================================
# parameter values to use in this script and the post-processing plotting script.
# =====================================

# ==========================================================================
# Parameters you might want to modify
# ==========================================================================
analytic_solution = 'raymond'  # can be 'raymond' or 'schoof'
kinflag = 1    # 1=apply kinematic bc (analytic soln) at points in the domain (discussed further below); 0=the run will be doubly periodic (preferred)
fillInitialGuess = 0  # 1=use the analytic solution as the initial guess for the velocity solver to speed convergence; 0=use the default 0-velocity initial guess


# Domain parameters
streamHalfWidth = 25000.0   # ice stream half-width, in m - used for both raymond & schoof formulations
alongFlowLength = 30000.0   # the desired along-flow length of the domain, in m; set to -1 to get a square domain
H = 1000.0       # ice thickness
dsdx = -1.0e-3   # bed (and surface) slope in the x-direction (y-direction is flat)

# Physical parameters
rho = 910.0   # ice density kg/m3
g = -9.81     # gravity m/s2
n = 3         # flow law exponent
A = 1e-16     # flow rate factor in Pa^-3 yr^-1
# ==========================================================================


# ======================================
# -- Functions for analytic solutions --
# ======================================

# raymond yield stress
def raymond_tau(yy):
  tau0 = 5.2e3*numpy.ones(yy.shape)         # set the stream value everywhere
  tau0[numpy.absolute(yy)>=streamHalfWidth] = 0.7e5        # set a very large value  outside the stream
  return tau0

# raymond velocity solution
def raymond_uvel(yy):
  tau0r = raymond_tau(yy)
  tau0r[tau0r>taud] = taud
  ur = 2.0 * A / (n+1.0) * ( (taud - tau0r)/H )**n * ( streamHalfWidth**(n+1) - numpy.absolute(yy)**(n+1) )
  ur[ur<0.0] = 0.0
  return ur

# schoof solution parameters
m = 1.55  # schoof exponent
L = streamHalfWidth / (m+1.0)**(1.0/m)  # This comes from the line above eq. 4.3 in schoof (2006)

# schoof yield stress distribution
def schoof_tau(yy):
  return taud * numpy.absolute( yy / L )**m

# schoof velocity solution
def schoof_uvel(yy):
  B = A**(-1.0/n)
  us = -2.0 * taud**3 * L**4 / (B**3 * H**3) * ( ((yy/L)**4 - (m+1.0)**(4.0/m))/4.0 - 3.0*( numpy.absolute(yy/L)**(m+4.0) \
    - (m+1.0)**(1.0+4.0/m) )/((m+1.0)*(m+4.0)) + 3.0*( numpy.absolute(yy/L)**(2.0*m+4.0) - (m+1.0)**(2.0+4.0/m) )/((m+1.0)**2*(2.0*m+4.0)) \
    - ( numpy.absolute(yy/L)**(3.0*m+4.0) - (m+1.0)**(3.0+4.0/m) )/ ( (m+1.0)**3*(3.0*m+4.0)) )

  # Some adjustments to the analytic profile - not entirely sure why these are needed.
  ind = numpy.nonzero( numpy.absolute(yy) >= streamHalfWidth )
  us[ind] = 0.0

  return us


# ======================================
# -- Values derived from things above --
# ======================================
taud = rho * g * H * dsdx  # Driving stress
# Calculate a good size for the size of the domain outside of the stream (in m)
if analytic_solution == 'raymond':
  strongWidth = 5.0 * H  # 5 ice thicknesses should get us beyond the zone of lateral stress transfer.  Adjust as needed
elif analytic_solution == 'schoof':
  # schoof (2006) uses a domain size that is 3L on either side of the central axis
  strongWidth = 3.0 * L - streamHalfWidth

# Calculating the actual domain size will happen later after the config file is read



# =====================================
# The actual script to run
# =====================================
if __name__ == '__main__':

    # =====================================
    # Initialization stuff

    # Parse options
    from optparse import OptionParser
    optparser = OptionParser()
    optparser.add_option("-c", "--config", dest="configfile", type='string', default='stream.config.in', help="Name of .config file to use to setup and run the test case", metavar="FILE")
    optparser.add_option('-m','--parallel',dest='parallel',type='int', help='Number of processors to run the model with: if specified then execute run in parallel', metavar="NUMPROCS")
    optparser.add_option('-e','--exec',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
    optparser.add_option('-s','--stream-size',dest='stream_grid_size',default=25,type='int',help='Number of cells to use to model the ice stream portion of the domain (values <19 may not work properly for all problems).')
    optparser.add_option('-v','--vert-grid-size',dest='vertical_grid_size',default=2,type='int',help='Number of vertical layers to use (upn); minimum value = 2')

    for option in optparser.option_list:
        if option.default != ("NO", "DEFAULT"):
            option.help += (" " if option.help else "") + "[default: %default]"
    options, args = optparser.parse_args()

    # Other modules needed by the main script
    from netCDF import *
    import ConfigParser, sys, os

    # =====================================
    # Create a netCDF file according to the information in the config file.
    try:
        configParser = ConfigParser.SafeConfigParser()
        configParser.read(options.configfile)
        nz = int(configParser.get('grid','upn'))
        filename = configParser.get('CF input', 'name')
    except:
        sys.exit('Error parsing ' + options.configfile)


    # =====================================
    # Figure out domain size information
    nStream = int(options.stream_grid_size)
    dy = 2.0 * streamHalfWidth / float(nStream)
    nStrongStrip = int(round(strongWidth / dy))  # Figure out the number of cells we need to add to get as close the the desired width of the strong region as possible (note: may want to use ceil() instead of round() here)
    ny = nStream + 2 * nStrongStrip  # a +1 is needed to convert from y0 to y1 but we leaving it off lets the stream boundaries fall on the y0 grid, which is needed to best match the analytic solution
    dx = dy  # always want this
    print 'Number of cells for stream (N-S):', nStream
    print 'Number of cells for entire domain (N-S):', ny
    print 'dy=dx=',dy
    print 'Domain N-S (across-flow) width (m):', ny*dy
    if alongFlowLength < 0:
      nx = ny  # square domain
    else:
      nx = int(round(alongFlowLength / dx))
    print 'Domain E-W (along-flow) width (m):', nx*dx

    # Now set the values in our internal config object to get written out later.
    configParser.set('grid','ewn',str(nx))
    configParser.set('grid','nsn',str(ny))
    configParser.set('grid','dew',str(dx))
    configParser.set('grid','dns',str(dy))

    if options.vertical_grid_size:
      nz = int(options.vertical_grid_size)
      configParser.set('grid','upn',str(nz))
      print 'Number of vertical layers has been set to:', nz

    # Check domain sizes for usefulness
    if (nStream % 2) == 0 and analytic_solution == 'schoof':
      print "Warning: For the schoof version, you might want the number of cells in the stream to be an odd number."

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

    x1 = dx*numpy.arange(nx,dtype='float64')
    y1 = dy*numpy.arange(ny,dtype='float64') - dy*float(ny-1)/2.0  # make the y-coordinates centered about 0
    x0 = dx/2.0 + x1[:-1] # staggered grid
    y0 = dy/2.0 + y1[:-1]

    # Make sure the edge of the stream lands on the grid cells on the y0 grid.  This should always happen with the logic above, so this check should never be activated.
    if (analytic_solution == 'raymond') and (not True in (numpy.absolute(streamHalfWidth-y0) < 0.0001)):
      sys.exit('Error: the stream edge does not land on the y0 grid so the stream will not be resolved adequately for the raymond case.  Adjust your domain size and/or stream size and/or horizontal resolution.  \nstream half width='+str(streamHalfWidth)+'\ny0 grid has values at:\n '+str(y0[:]))

    # Make sure we have at least two non-stream rows on each side
    if (numpy.absolute(y0[:])>streamHalfWidth).sum() < 4:
      sys.exit('Error: there are less than two non-stream rows on each side of the stream.  Adjust your domain size and/or stream size and/or horizontal resolution.  \nstream half width='+str(streamHalfWidth)+'\ny0 grid has cell centers at:\n '+str(y0[:]))

    # If all is good, continue.
    netCDFfile.createVariable('time','f',('time',))[:] = [0]
    netCDFfile.createVariable('x1','f',('x1',))[:] = numpy.float32(x1)
    netCDFfile.createVariable('y1','f',('y1',))[:] = numpy.float32(y1)
    netCDFfile.createVariable('x0','f',('x0',))[:] = numpy.float32(x0)
    netCDFfile.createVariable('y0','f',('y0',))[:] = numpy.float32(y0)

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
        tau0Profile = raymond_tau(y0)
        uvelProfile = raymond_uvel(y0)
    elif analytic_solution == 'schoof':
        tau0Profile = schoof_tau(y0)
        uvelProfile = schoof_uvel(y0)
    else:
        sys.exit("Error: Invalid value for 'analytic_solution'.")
    for i in range(nx-1):
      tauf[0,:,i] = tau0Profile


    # =======================================
    # Save the required variables to the netCDF file.
    netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg
    netCDFfile.createVariable('tauf','f',('time','y0','x0'))[:] = tauf

    if kinflag == 1 or fillInitialGuess == 1:
        netCDFfile.createVariable('uvel','f',('time','level','y0','x0'))
        netCDFfile.createVariable('vvel','f',('time','level','y0','x0'))

    if kinflag == 1:
        # setup Dirichlet boundary conditions for uvel and/or vvel at points in the domain

        dudy = numpy.gradient( uvelProfile, dy )
        vvelProfile = -dudy*dy

        kinbcmask = numpy.zeros([1,ny-1,nx-1],dtype='int32')
        uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
        vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')


    # =================================================================
    # fill both uvel and vvel at the upstrean and downstream domain ends

        # Fill first column
#        i = 0
#        uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
#        vvel[0,:,:,i] = -numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical
#        kinbcmask[0,:,i] = 1

        # Fill last column
#        i = nx-1 - 1
#        uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
#        vvel[0,:,:,i] = numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical
#        kinbcmask[0,:,i] = 1

    # =================================================================
    # fill both uvel and vvel at the upstrean and downstream domain ends
        # Fill just a single across-flow profile in domain interior
        i = 2
        uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
#        vvel[0,:,:,i] = -numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical
        kinbcmask[0,:,i] = 1

        netCDFfile.variables['uvel'][:] = uvel[:]
        netCDFfile.variables['vvel'][:] = vvel[:]
        netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kinbcmask[:]

    if fillInitialGuess == 1:
        # Fill the analytic solution into the initial guess to speed convergence
        dudy = numpy.gradient( uvelProfile, dy )
        vvelProfile = -dudy*dy

        uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
        vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

        for i in range(nx-1):
          uvel[0,:,:,i] = numpy.tile(uvelProfile, [nz, 1])  # uniform in the vertical
          vvel[0,:,:,i] = numpy.tile(vvelProfile, [nz, 1])  # uniform in the vertical

        netCDFfile.variables['uvel'][:] = uvel[:]
        netCDFfile.variables['vvel'][:] = vvel[:]

    netCDFfile.close()


    # =======================================
    # Update periodic offset in config file
    offset = -dsdx * dx * nx
    configParser.set('parameters', 'periodic_offset_ew', str(offset))

    # Save periodic offset and possibly nsn, ewn, dns, dew, upn to the config file
    runConfigFile = 'stream.config'
    configFile = open(runConfigFile,'w')
    configParser.write(configFile)
    configFile.close()



    # =====================================
    # Run CISM

    print 'Running CISM'
    print '============\n'
    if options.parallel == None:
       # Perform a serial run
       os.system(options.executable + ' ' + runConfigFile)
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
             sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run the model manually with something like: mpirun -np 4 ./cism_driver stream.config')
          runstring = mpiexec + ' ' + options.executable + ' ' + runConfigFile
          print 'Executing parallel run with:  ' + runstring + '\n\n'
          os.system(runstring)  # Here is where the parallel run is actually executed!


