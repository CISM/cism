#!/usr/bin/env python

# This script runs any or all of the various MISMIP3d experiments.
# See this paper for details about the MISMIP3d experiment:
# Pattyn et al., Grounding-line migration in plan-view marine ice-sheet models: results of the ice2sea MISMIP3d intercomparison, J. Glaciol., 59, doi:10.3189/2013JoG12J129, 2013.

import sys, os
import shutil
import fileinput
import numpy as np
from optparse import OptionParser
from ConfigParser import ConfigParser
from netCDF4 import Dataset


###################################
# string constants and parameters #
###################################

templateJobScriptFile = 'runCISM.cheyenne.template'  # template job submission script


####################################
# Function used later in the code #
####################################


# This function returns the basal friction perturbation for experiment P75S.
# Note: x, y and xGL (the grounding line position) need to be in km.
def computeBasalPerturbation(x,y,xGL):
    a   = 0.75  # perturbation amplitude
    xc  = 150   # in km
    yc  = 10    # in km
    yGL = 50    # in km (The function is center at the center line of the domain)

    # Calculating the spatial variation.
    xpart    = -(x-xGL)**2/(2*xc**2)
    ypart    = -(y-yGL)**2/(2*yc**2)
    Cspatial = 1-a*np.exp(xpart + ypart)
    return Cspatial


# This function is used to replace a string by another in a file.
def replace_string(file, oldString, newString):
    # Replace the first instance of oldString with newString in the input file.
    stringFound = False
    for line in fileinput.FileInput(file, inplace=True):
        if (not stringFound) and (oldString in line):
            line = line.replace(oldString, newString)
            stringFound = True
            print line,
        else:
            print line,



########
# Code #
########


# Parse options.
optparser = OptionParser()

optparser.add_option('-e', '--exec',     dest='executable', type = 'string', default='./cism_driver', help='Path to the CISM executable')
optparser.add_option('-x', '--expt',     dest='experiment', type='string',   default='all', help='MISMIP3d experiment(s) to run', metavar='EXPT')
optparser.add_option('-n', '--parallel', dest='parallel',   type='int', help='Number of processors: if specified then run in parallel', metavar='NUMPROCS')
optparser.add_option('--job',     action="store_true", dest='jobscript', help='option to set up a script to run on HPC')
optparser.add_option('--submit',  action="store_true", dest='jobsubmit', help='option to submit the run script directly after executing this script')


for option in optparser.option_list:
    if option.default != ('NO', 'DEFAULT'):
        option.help += (' ' if option.help else '') + '[default: %default]'
options, args = optparser.parse_args()

if options.experiment == 'all':
    experiments = ['Stnd', 'P75S', 'P75R']
    print 'Run all the MISMIP3d experiments'
elif options.experiment == 'allP75':
    experiments = ['P75S', 'P75R']
    print 'Run all the MISMIP3d P75 (S) and (R) experiments'
elif options.experiment in ['Stnd', 'P75S', 'P75R']:
    experiments = [options.experiment]
    print 'Run experiment', options.experiment
else:
    sys.exit('Please specify experiment(s) from this list: Stnd, P75S, P75R')


if (options.jobsubmit) and (not options.jobscript):
    print 'You are trying to submit a batch job without the job option. No job will be submitted.'

if (options.parallel is not None) and (options.jobscript):
    print 'The number of processors will not be taken into account as you will be submitting a batch job script.'

restart = 0  # constant meant to keep track of restart status.

# Loop through experiments.
for expt in experiments:

    # Change to directory for this experiment.
    os.chdir(expt)
    print 'Changed directory to ', expt

    # Set and open config file for reading.
    configfile = 'mismip3d' + expt + '.config'
    config = ConfigParser()
    config.read(configfile)

    inputFile   = config.get('CF input',   'name')
    outputFile  = config.get('CF output',  'name')
    outputFreq  = config.get('CF output',  'frequency')
    endTime     = config.get('time',       'tend')
    endTime     = float(endTime)
    
    # Time buffer to remedy the restart pb when switching to different time step within a run.
    buffer = float(outputFreq) - 1.

    # Read output file content information.
    lastTimeEntry     = 0
    lastEntryInternal = 0
    sizeTimeOutput    = 0
    if os.path.exists(outputFile):
        outputData = Dataset(outputFile,'r')
        if outputData['time'].size != 0:
            sizeTimeOutput = outputData['time'].size
            lastTimeEntry  = outputData['time'][-1]
        
        outputData.close()

    # Take action based on the stage of the experimental run.
    if (lastTimeEntry >= (endTime - buffer)) and (sizeTimeOutput > 1):
        # The run for this A value is done, moving to the next one.
        pass
    elif (lastTimeEntry < endTime) and (sizeTimeOutput > 1):
        # The run for this A value is not done and needs to continue.
        
        # Make sure restart is set to 1 in config file.
        config.set('options', 'restart', 1)
        
        # Write to config file.
        with open(configfile, 'w') as newconfigfile:
            config.write(newconfigfile)
    
    else:
        print 'there is nothing to restart from, executing from the beginning'


    # Make sure we are starting from the final time slice of the input file,
    # (Except for Stnd, the input file is a restart file from a previous run.)



    if expt != 'Stnd':

        # Edit the 'time' entry in the [CF input] section.
        inputfile = config.get('CF input', 'name')
        try:
           inputfileOut = config.get('CF input', 'nameout')
        except:
           pass  # no nameout name in config file
        
        infile = Dataset(inputfile,'r')
        ntime = len(infile.dimensions['time']) # reading the last time slice of the previous experiment
        config.set('CF input', 'time', ntime)
        print 'input file =', inputfile
        print 'time slice =', ntime
       
        infile.close()


        # For experiment P75S we need to induce the perturbation and modify C_space_factor.
        if (expt == 'P75S') and (restart==0):
            
            # We need to calculate the grounding line location from the input file of the experiment at the center line location.
            infileOut = Dataset(inputfileOut,'r')
            ntimeout  = len(infileOut.dimensions['time']) # reading last time slice of previous ext output file
            f_ground  = infileOut.variables['f_ground'][:,:,:]
            infileOut.close()

            print 'opening inputfile ',inputfile
            infile = Dataset(inputfile,'r+')
            nx = len(infile.dimensions['x0'])
            ny = len(infile.dimensions['y0'])
            x0 = infile.variables['x0'][:]
            y0 = infile.variables['y0'][:]

            ycenter = int(ny/2)
            print 'ycenter=',ycenter
                
            for i in range(1,nx-1):
                if (f_ground[ntimeout-1,ycenter,i] > 0) and (f_ground[ntimeout-1,ycenter,i+1] == 0):
                     xGL = x0[i]
                     print 'grounding line at ycenter =', xGL
        
            # We now write C_space factor to the file.
            print 'writing peturbation to C_space_factor'
            for i in range(nx):
                for j in range(ny):
                     infile.variables['C_space_factor'][-1,j,i] = computeBasalPerturbation(x0[i]/1000,y0[j]/1000,xGL/1000)

            # We first check that the variable "internal_time" from the Stnd restart file has 0 for its last entry.
            lastentry = infile['internal_time'][-1]
            if lastentry != 0:
                infile['internal_time'][:] = infile['internal_time'][:] - lastentry
                print 'the new internal_time array is ', infile['internal_time'][:]
            
            infile.close()

    
        # We need to reset C_space_factor to 1 and internal_time to 0.
        if (expt == 'P75R') and (restart==0):
             infile = Dataset(inputfile,'r+')
             print 'opening inputfile ',inputfile
             infile.variables['C_space_factor'][:,:,:] = 1

             lastentry = infile['internal_time'][-1]
             if lastentry != 0:
                 infile['internal_time'][:] = infile['internal_time'][:] - lastentry
                 print 'the new internal_time array is ', infile['internal_time'][:]

             infile.close()
        

        # Write the modified config file.
        with open(configfile, 'w') as newconfigfile:
            config.write(newconfigfile)



    # Checking whether we want to submit a job script on HPC.
    if options.jobscript:
        
        # New job submission script for current experiment.
        newJobScriptFile = 'runCISM.cheyenne.' + expt

        # Make a copy of the job submission script.
        try:
            shutil.copy('../' +  templateJobScriptFile, newJobScriptFile)
        except:
            sys.exit('Could not copy ' +  templateJobScriptFile)
    
        print 'Created master job submission script ', newJobScriptFile, ' for experiment ', expt


        # Modifying the walltime and computing power based on resolution (and Bill's experience).
        walltimeLook = '#PBS -l walltime=00:01:00'
        HPCpowerLook = '#PBS -l select=1:ncpus=1:mpiprocs=1'

        # Obtain the resolution from config file.
        res = config.get('grid', 'dew')

        if res=='8000.0':
            if (expt=='Stnd') or (expt=='P75R'):
                walltimeReplace = '#PBS -l walltime=00:05:00'
                HPCpowerReplace = '#PBS -l select=1:ncpus=36:mpiprocs=36'
            else:
                walltimeReplace = '#PBS -l walltime=00:01:00'
                HPCpowerReplace = '#PBS -l select=1:ncpus=36:mpiprocs=36'

        elif res=='4000.0':
            if (expt=='Stnd') or (expt=='P75R'):
                walltimeReplace = '#PBS -l walltime=00:30:00'
                HPCpowerReplace = '#PBS -l select=1:ncpus=36:mpiprocs=36'
            else:
                walltimeReplace = '#PBS -l walltime=00:01:00'
                HPCpowerReplace = '#PBS -l select=1:ncpus=36:mpiprocs=36'

        elif res=='2000.0':
            if (expt=='Stnd') or (expt=='P75R'):
                walltimeReplace = '#PBS -l walltime=02:00:00'
                HPCpowerReplace = '#PBS -l select=2:ncpus=36:mpiprocs=36'
            else:
                walltimeReplace = '#PBS -l walltime=00:08:00'
                HPCpowerReplace = '#PBS -l select=1:ncpus=36:mpiprocs=36'

        elif res=='1000.0':
            if (expt=='Stnd') or (expt=='P75R'):
                walltimeReplace = '#PBS -l walltime=06:00:00'
                HPCpowerReplace = '#PBS -l select=3:ncpus=36:mpiprocs=36'
            else:
                walltimeReplace = '#PBS -l walltime=00:05:00'
                HPCpowerReplace = '#PBS -l select=3:ncpus=36:mpiprocs=36'

        elif res=='500.0':
            if (expt=='Stnd') or (expt=='P75R'):
                walltimeReplace = '#PBS -l walltime=10:00:00'
                HPCpowerReplace = '#PBS -l select=4:ncpus=36:mpiprocs=36'
            else:
                walltimeReplace = '#PBS -l walltime=00:05:00'
                HPCpowerReplace = '#PBS -l select=4:ncpus=36:mpiprocs=36'

        elif res=='250.0':
            if (expt=='Stnd') or (expt=='P75R'):
                walltimeReplace = '#PBS -l walltime=10:00:00'
                HPCpowerReplace = '#PBS -l select=6:ncpus=36:mpiprocs=36'
            else:
                walltimeReplace = '#PBS -l walltime=00:05:00'
                HPCpowerReplace = '#PBS -l select=6:ncpus=36:mpiprocs=36'

        else:
            sys.exit('You are running with a customized resolution or one that is not supported. Adjust the job script manually')


        # Replacing computing resources to the job script file.
        try:
            linecheck = walltimeLook
            replace_string(newJobScriptFile, walltimeLook, walltimeReplace)
            linecheck = HPCpowerLook
            replace_string(newJobScriptFile, HPCpowerLook, HPCpowerReplace)
        except:
            print 'this line does not exist: ', linecheck
            
        # Modifying the executable line based on config file name.
        execlineLook    = 'mpiexec_mpt ./cism_driver mismip3d.config'
        execlineReplace = 'mpiexec_mpt ./cism_driver ' + configfile
        try:
            replace_string(newJobScriptFile, execlineLook, execlineReplace)
        except:
            print 'this line does not exist: ', execlineLook

        if options.jobsubmit:
            # Submitting the job.
            try:
                submitscriptstring = 'qsub' + ' ' + newJobScriptFile
                print 'submitting the batch job script '
                os.system(submitscriptstring)
            except:
                sys.exit('Could not submit job. Check if you are running this script on HPC Cheyenne.')

            sys.exit('Interrupting python script execution since job got submitted')

        else:
            sys.exit('Interrupting python script execution since experiment will be run by submitting a batch job script')



    # Run CISM (only if no job script is being created and/or run).

    print 'parallel =', options.parallel

    if options.parallel == None:
        # Perform a serial run.
        os.system(options.executable + ' ' + configfile)
    else:
        # Perform a parallel run.
        if options.parallel <= 0:
            sys.exit( 'Error: Number of processors specified for parallel run is <=0.' )
        else:
            # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure.
            if os.system('which openmpirun > /dev/null') == 0:
                mpiexec = 'openmpirun -np ' + str(options.parallel)
            elif os.system('which mpirun > /dev/null') == 0:
                mpiexec = 'mpirun -np ' + str(options.parallel)
            elif os.system('which aprun > /dev/null') == 0:
                mpiexec = 'aprun -n ' + str(options.parallel)
            elif os.system('which mpirun.lsf > /dev/null') == 0:
                # mpirun.lsf does NOT need the number of processors (options.parallel).
                mpiexec = 'mpirun.lsf'
            else:
                sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run manually with something like: mpirun -np 4 ./cism_driver mismip3dInit.config')

            runstring = mpiexec + ' ' + options.executable + ' ' + configfile
            print 'Executing parallel run with:  ' + runstring + '\n\n'

            # Here is where the parallel run is actually executed.
            os.system(runstring)

    print 'Finished experiment', expt

    # Change to parent directory and continue.
    os.chdir('..')
