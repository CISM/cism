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

# Raymond yield stress
def raymond_tau(yy):
  return 5.2e3*numpy.ones(yy.shape)

# Schoof yield stress distribution
def schoof_tau(yy):
  return taud * numpy.absolute( yy_centered / L )**m



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
    for option in optparser.option_list:
        if option.default != ("NO", "DEFAULT"):
            option.help += (" " if option.help else "") + "[default: %default]"
    options, args = optparser.parse_args()


    from netCDF import *
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

    if (ny % 2) != 0:
      sys.exit("Error: 'nsn' in the config file must be divisible by 2")

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


    y0_centered = y0 - ny * dy / 2.0  # use y coord that are symmetrical about 0

    if analytic_solution == 'raymond':
        tau0profile = raymond_tau(y0_centered)
    elif analytic_solution == 'schoof':
        tau0profile = schoof_tau(y0_centered)
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



##################################
# old stuff from the .m script that may be useful in the future

###if( flag == 0 )         % assign Raymond profile
###    tauf_profile = tau0r*ones(r-7,1);        
###    tauf = 0.7e5 * ones( r-1, c-1 );
###    if( kinflag ~= 1 )
###        tauf(4:end-3,3:end-2) = repmat( tauf_profile, 1, c-5 );     %% no slip at/near up/downstream ends
###    else
###        tauf(4:end-3,:) = repmat( tauf_profile, 1, c-1 );     %% use periodic bcs for cont. along flow
###    end
###    u_profile = repmat( [ 0 0 fliplr(ur) ur(2:end) 0 0 ], levels, 1 );     
###else                    % assign Schoof
###    
###    tauf_profile = [ fliplr(tau0s(2:end)), tau0s ]';        
###    tauf = 0.7e5 * ones( r-1, c-1 );
###    if( kinflag ~= 1 )
###        tauf(3:end-2,3:end-2) = repmat( tauf_profile, 1, c-5 );     %% no slip at/near up/downstream ends
###    else
###        tauf(3:end-2,:) = repmat( tauf_profile, 1, c-1 );     %% use periodic bcs for cont. along flow
###    end
###    u_profile = repmat( [ 0 0 fliplr(us) us(2:end) 0 0 ], levels, 1 );
###end
###%% use shear strain rate at up/downstream ends to calculate across-flow velocity profile
###%% necessary to balance the shear
###x = [0:dx:dx*length(u_profile)-1];
###dudy = gradient( u_profile(1,:), dy );
###v_profile = -dudy*dy;

###figure(999), clf
###subplot( 2,1,1 ),hold on
###plot( x, dudy, 'bo-' ), box on
###xlabel( 'dist (m)' ), ylabel( 'shear strain rate (1/a)' )
###subplot( 2,1,2 ),hold on
###plot( x, v_profile, 'bo-' ), box on
###xlabel( 'dist (m)' ), ylabel( 'across-flow component of vel (m/a)' )

###v_profile = repmat( v_profile, levels, 1 );

###% figure(200),clf
###% subplot(2,2,1),hold on
###% imagesc( thck ), axis xy, axis equal, axis tight, colorbar, title( 'thickness (m)' )
###% subplot(2,2,2),hold on
###% imagesc( usrf ), axis xy, axis equal, axis tight, colorbar, title( 'surface (m)' )
###% subplot(2,2,3),hold on
###% imagesc( topg ), axis xy, axis equal, axis tight, colorbar, title( 'bed (m)' )
###% subplot(2,2,4),hold on
###% imagesc( tauf/1e3 ), axis xy, axis equal, axis tight, colorbar, title( 'yield stress (kPa)' )

###%% optional: add analytic solution at up/downstream ends as kin vel bc
###kinbcmask = zeros(size(tauf));
###uvelhom = zeros( levels, r-1, c-1 );
###vvelhom = zeros( levels, r-1, c-1 );

###if( kinflag == 1)
###    uvelhom( :, :, end ) = u_profile; uvelhom( :, :, 1 ) = u_profile; 
###    vvelhom( :, :, end ) = v_profile; vvelhom( :, :, 1 ) = -v_profile; 
###    kinbcmask(:,1) = 1; kinbcmask(:,end) = 1;
###    ind = find( tauf <= 0 ); tauf(ind) = 1e-10;
###end

###%% for newer code, uvelhom = uvel, etc.
###uvel = uvelhom; vvel = vvelhom;

###save stream.mat usrf topg thck tauf levels 

###%% spit out some other vars needed in the .config file
###periodic_offset = dx*(c-0.5)*dsdx;
###dx
###periodic_offset
