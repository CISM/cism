#!/usr/bin/env python2

#FIXME: More detailed description of this test case!!!
"""
Run an experiment with an confined ice shelf. 
"""

# Authors
# -------
# Written by Glen Granzow at the University of Montana on April 8, 2010
# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 release
# Reconfigured by Joseph H Kennedy at ORNL on April 27, 2015 to work with the regression testing


import os
import sys
import errno
import subprocess
import ConfigParser 

import numpy
import netCDF
from math import sqrt, exp

# Parse the command line options
# ------------------------------
import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

parser.add_argument('-c','--config', default='./shelf-confined.config', 
        help="The configure file used to setup the test case and run CISM")
parser.add_argument('-e','--executable', default='./cism_driver', 
        help="The CISM driver")
parser.add_argument('--hpc', nargs='?', const='aprun',
        help=" ".join(["Flag to Shortcut parallel run command lookup for High Performance Computing Systems.", 
                       "If flag apears without an argument, it will set run command to `aprun`,", 
                       "otherwise it will use the argument given."]))
parser.add_argument('-m', '--modifier', metavar='MOD', default='',
        help="Add a modifier to file names. FILE.EX will become FILE.MOD.EX")
parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
parser.add_argument('-o', '--output-dir', default='./output',
        help="Write all created files here.")
parser.add_argument('-q', '--quiet', action='store_true',
        help="Run the CISM process quietly.")
parser.add_argument('-s','--setup-only', action='store_true',
        help="Set up the test, but don't actually run it.")


# Additional test specific options:
#parser.add_argument('--scale', type=unsigned_int, default=0, 
#        help="Scales the problem size by 2**SCALE. SCALE=0 creates a 31x31 grid, SCALE=1 " 
#            +"creates a 62x62 grid, and SCALE=2 creates a 124x124 grid.")


# Some useful functions
# ---------------------

# function to make a directory, and not worry if it exists.
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# prep the command line functions
def prep_commands(args, config_name):
    driver = os.path.abspath(args.executable)
   
    quiet_mod = ''
    if args.quiet:
        quiet_mod = ' > '+config_name+'.oe'

    commands = []
    mkdir_p(args.output_dir)
    commands.append("cd "+os.path.abspath(args.output_dir))
    
    if args.hpc and (args.parallel > 0):
        mpiexec = args.hpc+' -n ' + str(args.parallel)+" "
    elif (args.parallel > 0):
        # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
        if os.system('which openmpirun > /dev/null') == 0:
            mpiexec = 'openmpirun -np ' + str(args.parallel)+" "
        elif os.system('which mpirun > /dev/null') == 0:
            mpiexec = 'mpirun -np ' + str(args.parallel)+" "
        elif os.system('which aprun > /dev/null') == 0:
            mpiexec = 'aprun -n ' + str(args.parallel)+" "
        elif os.system('which mpirun.lsf > /dev/null') == 0:
            # mpirun.lsf does NOT need the number of processors
            mpiexec = 'mpirun.lsf '
        else:
            print("Unable to execute parallel run!")
            print("   Please edit the script to use your MPI run command, or run the model manually with")
            print("   something like: mpirun -np 4 ./cism_driver shelf-confined.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands


# the main script function
# ------------------------
def main():
    """
    Run the confined shelf test.
    """
    
    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------
    #scale_factor = 2 ** args.scale
    
    try:
        config_parser = ConfigParser.SafeConfigParser()
        config_parser.read( args.config )
        
        nz = int(config_parser.get('grid','upn'))
        
        nx = int(config_parser.get('grid','ewn'))#*scale_factor
        ny = int(config_parser.get('grid','nsn'))#*scale_factor
        dx = float(config_parser.get('grid','dew'))#/float(scale_factor)
        dy = float(config_parser.get('grid','dns'))#/float(scale_factor)
        
        file_name = config_parser.get('CF input', 'name')
        root, ext = os.path.splitext(file_name)

    except ConfigParser.Error as error:
        print("Error parsing " + args.config )
        print("   "), 
        print(error)
        sys.exit(1)
    
    res = str(nx).zfill(4)
    if args.parallel > 0:
        mod = args.modifier+'.'+res+'.p'+str(args.parallel).zfill(3)
    else:
        mod = args.modifier+'.'+res
    
    file_name = root+mod+ext
    config_name = root+mod+'.config'
    out_name = root+mod+'.out'+ext

    # create the new config file
    # --------------------------
    if not args.quiet: 
        print("\nCreating config file: "+config_name)
    
    config_parser.set('grid', 'ewn', str(nx))
    config_parser.set('grid', 'nsn', str(ny))
    config_parser.set('grid', 'dew', str(dx))
    config_parser.set('grid', 'dns', str(dy))

    config_parser.set('CF input', 'name', file_name)
    config_parser.set('CF output', 'name', out_name)
    config_parser.set('CF output', 'xtype', 'double')
    
    with open(config_name, 'wb') as config_file:
        config_parser.write(config_file)


    # create the input netCDF file
    # ----------------------------
    if not args.quiet: 
        print("\nCreating confined shelf netCDF file: "+file_name)
    try:
        nc_file = netCDF.NetCDFFile(file_name,'w',format='NETCDF3_CLASSIC')
    except TypeError:
        nc_file = netCDF.NetCDFFile(file_name,'w')

    nc_file.createDimension('time',1)
    nc_file.createDimension('x1',nx)
    nc_file.createDimension('y1',ny)
    nc_file.createDimension('level',nz)
    #nc_file.createDimension('staglevel',nz-1)
    #nc_file.createDimension('stagwbndlevel',nz+1)
    nc_file.createDimension('x0',nx-1) # staggered grid 
    nc_file.createDimension('y0',ny-1)

    x = dx*numpy.arange(nx,dtype='float32')
    y = dx*numpy.arange(ny,dtype='float32')

    nc_file.createVariable('time','f',('time',))[:] = [0]
    nc_file.createVariable('x1','f',('x1',))[:] = x
    nc_file.createVariable('y1','f',('y1',))[:] = y
    nc_file.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
    nc_file.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

    # Calculate values for the required variables.
    thk  = numpy.zeros([1,ny,nx],dtype='float32')
    topg = numpy.zeros([1,ny,nx],dtype='float32')
    beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')


    # *SFP* this has been changed so that the default value for 'flwa' is the same 
    # as in the EISMINT-shelf test documentation, tests 3 & 4, found at:
    # http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html

    ## Check to make sure that the flow law parameter in the config file is correct.
    #default_flwa = float(config_parser.get('parameters','default_flwa'))
    #if default_flwa != 4.6e-18:
    #  print 'WARNING: The parameter default_flwa in',args.config,'should be 4.6e-18'
    #  print '         Currently it is',default_flwa

    # Calculate values for the required variables.
    thk  = numpy.zeros([1,ny,nx],dtype='float32')
    topg = numpy.zeros([1,ny,nx],dtype='float32')
    beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')
    kbc  = numpy.zeros([1,ny-1,nx-1],dtype='int')
    acab = numpy.zeros([1,ny,nx],dtype='float32') # *sfp* added acab field for prog. runs 
    zero = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
    uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
    vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

    # *SFP* added topg var so that areas of no slip are consistent w/ grounded ice on bedrock
    topg[:] = -2000.0

    # *SFP* changed to be in line w/ EISMINT-shelf tests 3&4 

    # shelf bc applied at bottom (DEFAULT FOR TEST CASE - other options below for testing bcs)
    thk[0,4:-2,2:-2] = 500.     
    kbc[0,ny-4:,:]  = 1
    kbc[0,:,:3] = 1
    kbc[0,:,nx-4:] = 1
    topg[0,ny-4:,:]  = -440 
    topg[0,:,:4] = -440
    topg[0,:,nx-4:] = -440

    # shelf bc applied at top    
    #thk[0,2:-4,2:-2] = 500.     
    #kbc[0,:3,:]  = 1
    #kbc[0,:,:3] = 1
    #kbc[0,:,nx-4:] = 1
    #topg[0,:4,:]  = -440
    #topg[0,:,:4] = -440
    #topg[0,:,nx-4:] = -440

    # shelf bc applied at right     ! NOTE that shelf is wider slightly wider in ns than in ew direction  
    #thk[0,2:-2,2:-4] = 500.     
    #kbc[0,:,:3]  = 1
    #kbc[0,:3,:] = 1
    #kbc[0,ny-4:,:] = 1
    #topg[0,:,:4]  = -440
    #topg[0,:4,:] = -440
    #topg[0,ny-4:,:] = -440

    # shelf bc applied at left     ! NOTE that shelf is wider slightly wider in ns than in ew direction  
    #thk[0,2:-2,4:-2] = 500.     
    #kbc[0,:,nx-4:]  = 1
    #kbc[0,:3,:] = 1
    #kbc[0,ny-4:,:] = 1
    #topg[0,:,nx-4:]  = -440
    #topg[0,:4,:] = -440
    #topg[0,ny-4:,:] = -440


    beta[0,:,:] = 0 

    acab[:] = 0.25
    acab[0,ny-3:,:]  = 0    # zero out accum at edges to avoid buildup where u=0
    acab[0,:,:3] = 0
    acab[0,:,nx-3:] = 0

    # *SFP* calculate stream profile for upstream end
    for i in range(nx-2):
        x = float( i ) / (nx-2) - 0.5
        vvel[0,:,ny-4,i] = -1.5e3 * 1/(2*3.141592654*0.125) * numpy.exp( -x**2 / (2*0.125**2) )

    # Create the required variables in the netCDF file.
    nc_file.createVariable('thk',      'f',('time','y1','x1'))[:] = thk
    nc_file.createVariable('acab',     'f',('time','y1','x1'))[:] = acab
    nc_file.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kbc.astype(numpy.int32)
    nc_file.createVariable('topg',     'f',('time','y1','x1'))[:] = topg
    nc_file.createVariable('beta',     'f',('time','y0','x0'))[:] = beta
    nc_file.createVariable('uvel',  'f',('time','level','y0','x0'))[:] = zero

    # *sfp* first option below adds ice stream vel profile for kin bc at upstream end
    # *sfp* ... comment out for standard test case
    #nc_file.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = vvel
    nc_file.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = zero

    nc_file.close()
    mkdir_p(args.output_dir)
    subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)


    # Run CISM
    # --------
    command_list = prep_commands(args, config_name)
    commands_all = ["# SHELF-CONFINED"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet: 
            print("\nRunning CISM confined shelf test")
            print(  "================================\n")

        process = subprocess.check_call(str.join("; ",command_list), shell=True)
   
        try:
            subprocess.check_call("cd "+args.output_dir+"; "+result_mv, shell=True)
        except subprocess.CalledProcessError:
            pass 

        try:
            subprocess.check_call("cd "+args.output_dir+"; "+timing_mv, shell=True)
        except subprocess.CalledProcessError:
            pass

        if not args.quiet: 
            print("\nFinished running the CISM confined shelf test")
            print(  "=============================================\n")
    else:
        run_script = args.output_dir+os.sep+root+mod+".run" 
        
        with open(run_script,'w') as run_file:
            run_file.write('#!/usr/bin/env bash \n \n')
            for command in commands_all:
                run_file.write(command+" \n")

        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            print("\nFinished setting up the CISM confined shelf test")
            print(  "================================================")
            print(  "   To run the test, use: "+run_script)

# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())

