#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script runs ISMIP-HOM experiments using Glimmer-CISM.
# Output files are written in the "output" subdirectory.
# The script loops over experiments performing the following three steps:
# 1. Create two input files, a configuration file and a netCDF input file.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Write a standard ISMIP-HOM file from the data in the netCDF output file.
# When finished, any additional files written by Glimmer are moved to the "scratch" subdirectory.

# After running this script, run plotISMIPHOM.py to plot the results.
# See the accompanying README file for more information.
# To see all command line options run: python runISMIPHOM.py --help
# Written March 2, 2010 by Glen Granzow at the University of Montana.

# OptionParser callback that splits a comma-separated list.
def appendToList(option,str_opt,value,parser):
  listOfValues = getattr(parser.values,option.dest)
  if listOfValues == None: listOfValues = list()
  listOfValues += value.split(',')
  setattr(parser.values,option.dest,listOfValues)

defaultExperiments = ['a','b']      # ['a','b','c','d']
defaultSizes = ['5','10','20','40','80','160']

if __name__ == '__main__':
  import os
  import glob
  import shutil
  import ConfigParser
  from optparse import OptionParser
  from math import tan, sin, pi, exp
  from netCDF import *

# Parse the command line arguments
  parser = OptionParser()
  parser.add_option('-e','--exp',dest='experiments',type='string',action='callback',callback=appendToList,help='Specify ISMIP-HOM experiments to run')
  parser.add_option('-s','--size',dest='sizes',type='string',action='callback',callback=appendToList,help='Specify domain sizes to run')
  parser.add_option('-g','--grid-size',dest='horizontal_grid_size',type='int',help='overrides ewn and nsn in config file')
  parser.add_option('-v','--vert-grid-size',dest='vertical_grid_size',type='int',help='overrides upn in config file')
  parser.add_option('-r','--run',dest='executable',default='./cism_driver',help='Set path to the CISM executable')
  parser.add_option('-p','--prefix',dest='prefix',default='cis1',help='Prefix to use for model output files')
  parser.add_option('-f','--format-only',dest='format_only',action='store_true',help='Generate the config and NetCDF input files only')
  parser.add_option('-m','--parallel',dest='parallel',type='int',help='if specified then execute run in parallel')
  parser.add_option('-c','--cyclic',dest='cyclic',action='store_true',default=False,help='if specified then all fields, including scalars, are truly periodic across the domain (NOT true for ismip-hom)')
  for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
  options, args = parser.parse_args()
# If the user didn't specify a list of experiments or domain sizes, run the whole suite
  if options.experiments == None: options.experiments = defaultExperiments
  if options.sizes == None: options.sizes = defaultSizes

# Loop over the experiments requested on the command line
  for experiment in options.experiments:

#   Loop over the sizes requested on the command line
    for size in map(int,options.sizes):

      if experiment == 'f': 
        if size != 100 or len(options.sizes) > 1:
          print 'NOTE: Experiment f uses a domain size of 100 km only'
          size = 100

#     Create a configuration file by copying an existing example
      configParser = ConfigParser.SafeConfigParser()
      configParser.read('ishom.config')
#     Set or get the number of horizontal grid points in each direction
      if options.horizontal_grid_size == None:
        nx = int(configParser.get('grid','ewn'))
        ny = int(configParser.get('grid','nsn'))
      else:
        nx = ny = options.horizontal_grid_size
        configParser.set('grid','ewn',str(nx))
        configParser.set('grid','nsn',str(ny))
#     Set the grid spacing in the horizontal directions
      dx = float(size)*1000.0/float(nx)
      dy = float(size)*1000.0/float(ny)
      configParser.set('grid','dew',str(dx))
      configParser.set('grid','dns',str(dy))
#     Specify the netCDF input and output filenames
#     All files will be written in the "output" subdirectory
      filename = os.path.join('output','ishom.'+experiment+'.'+str(size)+'km')
      configParser.set('CF input', 'name',filename+'.nc')
      configParser.set('CF output','name',filename+'.out.nc')
      configParser.set('CF default','title', 'ISMIP-HOM Experiment ' + experiment.capitalize() )
#     add the vertical offset specification to the config if running actual ismip-hom test cases
      if not options.cyclic:
        if experiment in ('a','b'):
          offset = float(size)*1000.0 * tan(0.5 * pi/180.0)
        elif experiment in ('c','d'):
          offset = float(size)*1000.0 * tan(0.1 * pi/180.0)
        elif experiment in ('f'):
          offset = float(size)*1000.0 * tan(3.0 * pi/180.0)
        configParser.set('parameters', 'periodic_offset_ew', str(offset))

      if experiment in ('c','d'):
        # These tests have beta passed in from the input file, so change option accordingly.
        configParser.set('ho_options', 'which_ho_babc', '5')

##      #Optional: if doing experiment C, one can alternatively use the ho_babc option setup for this test case rather than passing in a beta
##      if experiment in ('c'):
##        configParser.set('ho_options', 'which_ho_babc', '8')

      # For test case F we need to make a few additional adjustments to the config
      if experiment in ('f'):
        configParser.set('ho_options', 'which_ho_efvs', '0')  # Set to efvs to be the constant value of 2336041.42829 hardcoded for this option - this corresponds to the value needed for this test case
        configParser.set('time', 'dt', '2.2')  # 2.4 yr is the longest dt ok for diffusive CFL, when using dx=dy=2500.0
        # Need to run to steady-state...  
        # It's close to SS by 400 years, but there are some long-period oscillations that still appear out to 1000 yrs.  Not sure yet how much longer than that to eliminate those.
        configParser.set('time', 'tend', '400.0')
        configParser.set('CF output', 'variables', 'uvel vvel uvel_icegrid vvel_icegrid topg thk usurf wvel_ho velnorm efvs adv_cfl_dt diff_cfl_dt')  # Include flwa, efvs and the CFL variables to the output file
        configParser.set('CF output', 'frequency', '25.0')  # we don't want to output a whole lot of time levels, but want to be able to see we've reached SS.

#     Make additional changes if requested on the command line
      if options.vertical_grid_size != None:
        configParser.set('grid','upn',str(options.vertical_grid_size))
#     Write the new configuration file
      configFile = open(filename+'.config','w')
      configParser.write(configFile)
      configFile.close()

#     Create the netCDF input file needed by CISM
      if netCDF_module == 'netCDF4':
        netCDFfile = NetCDFFile(filename+'.nc','w',format='NETCDF3_CLASSIC')
      else:
        netCDFfile = NetCDFFile(filename+'.nc','w')
      netCDFfile.createDimension('time',1)
      netCDFfile.createDimension('x1',nx)   # unstaggered grid
      netCDFfile.createDimension('y1',ny)
      netCDFfile.createDimension('x0',nx-1) # staggered grid 
      netCDFfile.createDimension('y0',ny-1)
      time = netCDFfile.createVariable('time','f',('time',))
      x1   = netCDFfile.createVariable('x1','f',('x1',))
      y1   = netCDFfile.createVariable('y1','f',('y1',))
      x0   = netCDFfile.createVariable('x0','f',('x0',))
      y0   = netCDFfile.createVariable('y0','f',('y0',))
      thk  = netCDFfile.createVariable('thk' ,'f',('time','y1','x1'))
      topg = netCDFfile.createVariable('topg','f',('time','y1','x1'))
      if experiment in ('c','d'):
        unstagbeta = netCDFfile.createVariable('unstagbeta','f',('time','y1','x1'))
      time[0] = 0
      x1[:] = [(i+0.5)*dx for i in range(nx)] # unstaggered grid
      y1[:] = [(j+0.5)*dy for j in range(ny)]
      x0[:] = [(i+1)*dx for i in range(nx-1)] # staggered grid 
      y0[:] = [(j+1)*dy for j in range(ny-1)]

#     Generate the ice thickness, bed topography, and (sometimes) 
#     basal friction coefficient for the experiment

      thickness  = list()
      topography = list()
      basalFriction = list()

      xx = [(i+0.5)*dx for i in range(nx)]
      yy = [(j+0.5)*dy for j in range(ny)]

      if experiment in ('a','b'):
        if not options.cyclic:              
          alpha = 0.5 * pi/180
        else: # optional flags to allow for truly periodic domain setup
          alpha = 0.
        zz = [4000-x1[i]*tan(alpha) for i in range(nx)]
      elif experiment in ('c','d'):
#        if not options.cyclic:              
        alpha = 0.1 * pi/180
#        else:  # optional flags to allow for truly periodic domain setup
#          alpha = 0.
        zz = [1000-x1[i]*tan(alpha) for i in range(nx)]
      elif experiment == 'f':
        alpha = 3.0 * pi/180
        zz = [7000-x1[i]*tan(alpha) for i in range(nx)]
        xc = (xx[0]+xx[-1])/2
        yc = (yy[0]+yy[-1])/2
        a0 = 100
        sigma2 = 10000**2

      omega = 2*pi / (size*1000)      
      for y in yy:
        row = list()
        for x in xx:
          if experiment == 'a':
            row.append(1000 - 500*sin(omega*x)*sin(omega*y))
          elif experiment == 'b':
            row.append(1000 - 500*sin(omega*x))
          elif experiment == 'c':
#            if options.cyclic:    # for test case w/ truly periodic domain, add some non-zero topog to force flow 
#              row.append(1000 - 500*sin(omega*x)*sin(omega*y))
#            elif not options.cyclic:
            row.append(1000 + 1000*sin(omega*x)*sin(omega*y))
          elif experiment == 'd':
            row.append(1000 + 1000*sin(omega*x))
          elif experiment == 'f':
            row.append(1000 - a0*exp(-((x-xc)**2+(y-yc)**2)/sigma2))
        if experiment in ('a','b','f'):
          thickness.append(row)
          if not options.cyclic:        
            topography.append([z-t for (z,t) in zip(zz,row)])
          else:  # options to allow for truly periodic domain setup
            topography.append([z for z in zz])
        else:
#          basalFriction.append(row[:-1])
          basalFriction.append(row[:])
 
      if experiment in ('a','b','f'):
        thk [:] = thickness
        topg[:] = topography
      elif experiment in ('c','d'):
        thk [:] = ny*[nx*[1000]]
        topg[:] = ny*[zz]
#        if not options.cyclic:         # options to allow for truly periodic domain setup
        unstagbeta[:] = basalFriction[:]
      netCDFfile.close()

      if not options.format_only:

#       Run CISM (NOTE two options here.)
        print 'Running',options.executable,'for experiment',experiment.upper(),'with domain size',size,'km'
        if options.parallel != None:
           if options.parallel > 0:
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
                 sys.exit('Unable to execute parallel run.  Please edit the script to use your MPI run command.')
              # Now run the model
              runstring = mpiexec + ' ' + options.executable + ' ' + filename + '.config'
              print 'Executing parallel run with:  ' + runstring + '\n\n'
              exitCode = os.system(runstring)
           else:
              print 'Number of processors specified for parallel run is <=0.  Skipping the running of the model.'
              exitCode = 0
        else:
           # Perform a serial run
           exitCode = os.system(options.executable + ' ' + filename + '.config')

        if exitCode != 0:
          print 'Error: The CISM run for experiment '+experiment+' at size '+str(size)+' did NOT complete successfully!'

#     Experiment f should be run for one size (100 km) only
      if experiment == 'f': break

# Clean up by moving extra files written by CISM to the "scratch" subdirectory
# Look for files with extension "txt", "log", or "nc"
  for files in glob.glob('*.txt')+glob.glob('*.log')+glob.glob('*.nc'):
#   Delete any files already in scratch with these filenames 
    if files in os.listdir('scratch'):
      os.remove(os.path.join('scratch',files))
#   Move the new files to scratch
    shutil.move(files,'scratch')

  if options.cyclic:
     print '\nYou have specified the --cyclic flag which enforces true periodicity rather than the periodic slope of the ISMIP-HOM test specifications.'
     print '  Note: This flag is only applied to experiments A and B.'
     print '  Note: This model setup should NOT be compared to ISMIP-HOM results!'



