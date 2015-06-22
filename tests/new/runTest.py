#!/usr/bin/env python

#FIXME: Write a detailed description of this test case!!! This will
#       be display as the help text header.
"""
A generic test template. 
"""

# Authors
# -------
# Template by Joseph H Kennedy at ORNL on April 27, 2015

import os
import sys
import errno
import subprocess
import ConfigParser 

import numpy
import netCDF
from math import sqrt

# Parse the command line options
# ------------------------------
import argparse
parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#NOTE: This uses the doc string at the top of this file as the help description and
#      automatically generates a usage stament, showing the defaults for the 
#      arguments created below.

# small helper function so argparse will understand unsigned integers
def unsigned_int(x):
    x = int(x)
    if x < 1:
        raise argparse.ArgumentTypeError("This argument is an unsigned int type! Should be an integer greater than zero.")
    return x

#FIXME: correct the default config file name.
parser.add_argument('-c','--config', default='./TEST.config', 
        help="The configure file used to setup the test case and run CISM")
parser.add_argument('-e','--executable', default='./cism_driver', 
        help="The CISM driver")
parser.add_argument('--hpc', action='store_true',
        help="Shortcuts parallel run command lookup for High Performance Computing Systems. Will set run command to `time apirun -n N`.")
parser.add_argument('-m', '--modifier', metavar='MOD', default='',
        help="Add a modifier to file names. FILE.EX will become FILE.MOD.EX")
parser.add_argument('-n','--parallel', metavar='N', type=unsigned_int, default=0, 
        help="Run in parallel using N processors.")
parser.add_argument('-o', '--output-dir', default='./output',
        help="Write all created files here.")
parser.add_argument('-q', '--quiet', action='store_true',
        help="Run the cism process quietly.")
parser.add_argument('-s','--setup-only', action='store_true',
        help="Set up the test, but don't actully run it.")

#FIXME: Put any additional arguments needed here:
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
        mpiexec = 'aprun -n ' + str(args.parallel)+" "
    elif (args.parallel > 0):
        # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure
        if os.system('which openmpirun > /dev/null') == 0:
            mpiexec = 'openmpirun -np ' + str(args.parallel)+" "
        elif os.system('which mpirun > /dev/null') == 0:
            mpiexec = 'mpirun -np ' + str(args.parallel)+" "
        elif os.system('which aprun > /dev/null') == 0:
            mpiexec = 'aprun -n ' + str(args.parallel)+" "
        elif os.system('which mpirun.lsf > /dev/null') == 0:
            # mpirun.lsf does NOT need the number of processors (options.parallel)
            mpiexec = 'mpirun.lsf '
        else:
            print("Unable to execute parallel run!")
            print("   Please edit the script to use your MPI run command, or run the model mannually with")
            #FIXME: correct the default config file name.
            print("   something like: mpirun -np 4 ./cism_driver TEST.config")
            sys.exit(1)
    else:
        mpiexec = ''

    commands.append(mpiexec+driver+" "+config_name+quiet_mod)

    return commands

# the main script function
# ------------------------
def main():
    """
    Run the test.
    """
    
    # check that file name modifier, if it exists, starts with a '-'
    if not (args.modifier == '') and not args.modifier.startswith('-') :
        args.modifier = '-'+args.modifier
         
    # get the configuration
    # ---------------------
    #FIXME: If you need to parse the configurate file, do it here. Uncomment what you need
    #       below and make sure to add any config options you will need. 
    #
    #NOTE: Use this with the scale option above if needed. You can write your own logit here.
    #      scale_factor = 2 ** args.scale
    #
    #try:
    #    config_parser = ConfigParser.SafeConfigParser()
    #    config_parser.read( args.config )
    #    
    #    nz = int(config_parser.get('grid','upn'))
    #    
    #    nx = int(config_parser.get('grid','ewn'))*scale_factor
    #    ny = int(config_parser.get('grid','nsn'))*scale_factor
    #    dx = float(config_parser.get('grid','dew'))/float(scale_factor)
    #    dy = float(config_parser.get('grid','dns'))/float(scale_factor)
    #    
    #    file_name = config_parser.get('CF input', 'name')
    #    root, ext = os.path.splitext(file_name)
    #
    #except ConfigParser.Error as error:
    #    print("Error parsing " + args.config )
    #    print("   "), 
    #    print(error)
    #    sys.exit(1)
    #
    #NOTE: when using a scale factor, you'll likely want to modify the output files
    #      with the resolution. 
    #res = str(nx).zfill(4)
    #mod = args.modifier+'.'+res
    #
    #file_name = root+mod+ext
    #config_name = root+mod+'.config'
    #out_name = root+mod+'.out'+ext

    # create the new config file
    # --------------------------
    #FIXME: You should ALWAYS write a copy of the config file you are using. Uncomment
    #       what you need here. 
    #if not args.quiet: 
    #    print("Creating config file: "+config_name)
    #
    #config_parser.set('grid', 'ewn', str(nx))
    #config_parser.set('grid', 'nsn', str(ny))
    #config_parser.set('grid', 'dew', str(dx))
    #config_parser.set('grid', 'dns', str(dy))

    #config_parser.set('CF input', 'name', file_name)
    #config_parser.set('CF output', 'name', out_name)
    #config_parser.set('CF output', 'xtype', 'double')
    #
    #with open(config_name, 'wb') as config_file:
    #    config_parser.write(config_file)



    # create the input netCDF file
    # ----------------------------
    #FIXME: If you need to create an input netCDF file, do so here. Uncomment what you need
    #       below and make sure you parsed any variabled you will need out of the config file!
    #if not args.quiet: 
    #    print("Creating TEST netCDF file: "+file_name)
    #try:
    #    nc_file = netCDF.NetCDFFile(file_name,'w',format='NETCDF3_CLASSIC')
    #except TypeError:
    #    nc_file = netCDF.NetCDFFile(file_name,'w')

    #nc_file.createDimension('time',1)
    #nc_file.createDimension('x1',nx)
    #nc_file.createDimension('y1',ny)
    #nc_file.createDimension('level',nz)
    #nc_file.createDimension('staglevel',nz-1)
    #nc_file.createDimension('stagwbndlevel',nz+1)
    #nc_file.createDimension('x0',nx-1) # staggered grid 
    #nc_file.createDimension('y0',ny-1)

    #x = dx*numpy.arange(nx,dtype='float32')
    #y = dx*numpy.arange(ny,dtype='float32')

    #nc_file.createVariable('time','f',('time',))[:] = [0]
    #nc_file.createVariable('x1','f',('x1',))[:] = x
    #nc_file.createVariable('y1','f',('y1',))[:] = y
    #nc_file.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
    #nc_file.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]

    ## Calculate values for the required variables.
    #thk  = numpy.zeros([1,ny,nx],dtype='float32')
    #topg = numpy.zeros([1,ny,nx],dtype='float32')
    #artm = numpy.zeros([1,ny,nx],dtype='float32')

    #FIXME: Calculate the thickness (thk) here:
    #thk[:] = 100.

    #FIXME: Calculate the bed topography (topg) here:
    #topg[:] = 0.

    #FIXME: specify a sfc temperature field so that temperature evol. can be calc. if desired
    #artm[:] = -15.0

    #FIXME: Create the required variables in the netCDF file.
    #nc_file.createVariable('thk', 'f',('time','y1','x1'))[:] = thk
    #nc_file.createVariable('topg','f',('time','y1','x1'))[:] = topg
    #nc_file.createVariable('artm','f',('time','y1','x1'))[:] = artm 

    #FIXME: Calculate optional fields that could be added to the initial condition file.  
    #tempstag = numpy.zeros([1,nz+1,ny,nx],dtype='float32')
    #beta = numpy.zeros([1,ny-1,nx-1],dtype='float32')
    #nc_file.createVariable('tempstag','f',('time','stagwbndlevel','y1','x1'))[:] = tempstag 
    #nc_file.createVariable('beta','f',('time','y0','x0'))[:] = beta

    # close the new netCDF and config file, and move it to the output directory (if given)
    #nc_file.close()
    #mkdir_p(args.output_dir)
    #subprocess.check_call("cp *rilinosOptions.xml "+args.output_dir, shell=True)
    #subprocess.check_call("mv "+file_name+" "+args.output_dir, shell=True)
    #subprocess.check_call("mv "+config_name+" "+args.output_dir, shell=True)

    # Run CISM
    # --------
    #FIXME: correct the name of your test.
    command_list = prep_commands(args, config_name)
    commands_all = ["# TEST"+mod+" test"]
    commands_all.extend( command_list )
    
    result_mv = "mv results "+root+mod+".results 2>/dev/null"
    timing_mv = "for file in cism_timing*; do mv $file "+root+mod+".$file 2>/dev/null; done"
    commands_all.append(result_mv)
    commands_all.append(timing_mv)
    commands_all.append(" ")
    
    if not args.setup_only:
        if not args.quiet:
            #FIXME: correct the name of your test.
            print("Running CISM TEST test")
            print("======================\n")


        command_list = prep_commands(args, config_name)
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
            #FIXME: correct the name of your test.
            print("\nFinished running the CISM TEST test")
            print(  "===================================\n")
    else:
        run_script = args.output_dir+os.sep+root+mod+".run" 
        
        run_file = open(run_script,'w') 
        
        run_file.write('#!/bin/bash \n')
        for command in commands_all:
            run_file.write(command+" \n")

        run_file.close()
        os.chmod(run_script, 0o755)   # uses an octal number!

        if not args.quiet:
            #FIXME: correct the name of your test.
            print("\nFinished setting up the CISM TEST test")
            print(  "======================================")
            print(  "   To run the test, use: "+run_script)

   


# Run only if this is being run as a script.
if __name__=='__main__':
    
    # get the command line arguments
    args = parser.parse_args()
    
    # run the script
    sys.exit(main())
