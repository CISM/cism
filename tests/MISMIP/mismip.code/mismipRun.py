#!/usr/bin/env python

# This script runs the MISMIP (1d) experiment.
# See this paper for details:
# Pattyn et al., Results of the Marine Ice Sheet Model Intercomparison Project, MISMIP,
#    The Cryosphere, 6, 573-588, 2012, doi:10.5194/tc-6-573-2012.


import sys, os
import shutil
import fileinput
import numpy as np

from optparse import OptionParser
from netCDF4 import Dataset
from ConfigParser import ConfigParser



###############################
# Constants used in this code #
###############################


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





###############################
# Functions used in this code #
###############################


def launchCism(executable, configfile, parallel):
    # Run CISM (only if no job script is being created and/or run).
    print 'parallel =', parallel
    
    if parallel == None:
        # Perform a serial run.
        os.system(executable + ' ' + configfile)
    else:
        # Perform a parallel run.
        if parallel <= 0:
            sys.exit( 'Error: Number of processors specified for parallel run is <=0.' )
        else:
            # These calls to os.system will return the exit status: 0 for success (the command exists), some other integer for failure.
            if os.system('which openmpirun > /dev/null') == 0:
                mpiexec = 'openmpirun -np ' + str(parallel)
            elif os.system('which mpirun > /dev/null') == 0:
                mpiexec = 'mpirun -np ' + str(parallel)
            elif os.system('which aprun > /dev/null') == 0:
                mpiexec = 'aprun -n ' + str(parallel)
            elif os.system('which mpirun.lsf > /dev/null') == 0:
                # mpirun.lsf does NOT need the number of processors (options.parallel).
                mpiexec = 'mpirun.lsf'
            else:
                sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command, or run manually with something like: mpirun -np 4 ./cism_driver mismip3dInit.config')
        
            runstring = mpiexec + ' ' + executable + ' ' + configfile
            print 'Executing parallel run with:  ' + runstring + '\n\n'
            
            # Here is where the parallel run is actually executed!
            os.system(runstring)




########
# code #
########

# Parse options.
optparser = OptionParser()

optparser.add_option('-e', '--exec',    dest='executable',type = 'string',default='./cism_driver', help='Path to the CISM executable')
optparser.add_option('-n', '--parallel',dest='parallel',  type='int',help='Number of processors: if specified then run in parallel', metavar="NUMPROCS")
optparser.add_option('-x', '--expt',    dest='experiment',type='string',default='all',     help='MISMIP experiment set to run', metavar="EXPT")
optparser.add_option('-s', '--stat',    dest='StatChoice',type='string',default='advance', help='Experiment status set to run', metavar='STATUS')
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
    sys.exit('Please specify bed type from this list: linear, poly.')


if options.experiment == 'all':
    experiments = As
    Astat       = Astatus
    print 'Running all the MISMIP experiments'
elif options.experiment == 'advance':
    experiments = AsAdvance
    Astat       = AstatusAdvance
    print 'Running advance experiments'
elif options.experiment == 'retreat':
    experiments = AsRetreat
    Astat       = AstatusRetreat
    print 'Running retreat experiments'
elif options.experiment in As:
    # In this case there might be 2 possibilities, advance or retreat.
    experiments = [options.experiment]
    if options.StatChoice == 'retreat':
        Astat = ['retreat']
    else:
        Astat = ['advance']
    
    print 'Running experiment ', options.experiment
else:
    sys.exit('Please specify experiment(s) from this list: all, advance, retreat or a single value from Pattyn et al.2012.')



# loop through A values.
bedType   = options.bedtopo
countStat = -1                # counter to access status matrix. Taking into account zero-arrays indexing

for expt in experiments:
    countStat = countStat+1
    stat      = Astat[countStat]
    
    # Change to bed type directory.
    os.chdir(bedType)
    
    # Change to advance or retreat directory.
    os.chdir(stat)

    # Change to A value directory.
    os.chdir(expt)
    print 'changed directory to ', stat+'/'+expt

    # Name of the restart pointer file.
    restartPointer = 'mismip_' + expt + '.pointer'

    # Read information from config file.
    configfile = 'mismip_' + expt + '.config'
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
        print 'Continuing experiment from restart.'

        # Make sure restart is set to 1 in config file.
        config.set('options', 'restart', 1)

        # Write to config file.
        with open(configfile, 'w') as newconfigfile:
            config.write(newconfigfile)

        # launch CISM to restart the run.
        launchCism(options.executable, configfile, options.parallel)

        # Re-opening the output file and check if experiment was manually interrupted or not.
        outputData    = Dataset(outputFile,'r')
        lastTimeEntry = outputData['time'][-1]
        outputData.close()
        if (lastTimeEntry >= (endTime - buffer)):
            print 'Finished experiment', expt
        else:
            print 'Experiment interrupted.'
            sys.exit('Terminating the run.')
                
    else:
        # Start the experiment from time = 0.
        if (expt == As[0]) and (stat=='advance'):
            print 'First A-value, beginning from initial setup.'
        else:
            print 'Restarting from previous A-value restart file.'

            inputData          = Dataset(inputFile,'r+')
            lastEntryInternal  = inputData['internal_time'][-1]
            inputslice = 1
            if lastentry != 0:
                inputData['internal_time'][:] = inputData['internal_time'][:] - lastEntryInternal
                print 'the new internal_time array is ', inputData['internal_time'][:]

            inputData.close()
                 
            # Set config file.
            config.set('CF input', 'time', inputslice)


        # launch CISM.
        launchCism(options.executable, configfile, options.parallel)
                
        # Re-open the output file and check if experiment was manually interrupted or not.
        outputData    = Dataset(outputFile,'r')
        lastTimeEntry = outputData['time'][-1]
        outputData.close()
        if (lastTimeEntry >= (endTime - buffer)):
            print 'Finished experiment', expt
        else:
            print 'Experiment interrupted.'
            sys.exit('Terminating the run.')


    print('Switching back to original directory.')
    # Change to parent directory and continue.
    os.chdir('../../..')
