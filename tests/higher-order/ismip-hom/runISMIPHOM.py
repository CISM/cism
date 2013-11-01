#)!/usr/bin/env python
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
defaultSizes = ['160','80','40','20','10','5']

if __name__ == '__main__':
  import os
  import glob
  import shutil
  import ConfigParser
  from optparse import OptionParser
  from math import tan, sin, pi, exp
  from netCDF import *
  import numpy as np

# Parse the command line arguments
  parser = OptionParser()
  parser.add_option('-e','--exp',dest='experiments',type='string',action='callback',callback=appendToList,help='Specify ISMIP-HOM experiments to run')
  parser.add_option('-s','--size',dest='sizes',type='string',action='callback',callback=appendToList,help='Specify domain sizes to run')
  parser.add_option('-g','--grid-size',dest='horizontal_grid_size',type='int',help='(overrides ewn and nsn in config file)')
  parser.add_option('-v','--vert-grid-size',dest='vertical_grid_size',type='int',help='(overrides upn in config file)')
  parser.add_option('-r','--run',dest='executable',default='./simple_glide',help='Set path to the CISM executable (defaults to simple_glide)')
  parser.add_option('-p','--prefix',dest='prefix',default='cis1',help='Prefix to use for model output files (defaults to cis1)')
  parser.add_option('-f','--format-only',dest='format_only',action='store_true',help='Generate the config and NetCDF input files only')
  parser.add_option('-m','--parallel',dest='parallel',type='int',help='if specified then execute run in parallel')
  parser.add_option('-c','--cyclic',dest='cyclic',action='store_true',default=False,help='if specified then all fields, including scalars, are truly periodic across the domain (NOT true for ismip-hom)')
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
        configParser.set('parameters', 'periodic_offset_ew', str(offset))
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
        #beta = netCDFfile.createVariable('beta','f',('time','y0','x0'))
        beta = netCDFfile.createVariable('beta','f',('time','y1','x1'))
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
        zz = [6000-x1[i]*tan(alpha) for i in range(nx)]
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
#        beta[:] = basalFriction[:-1]
        beta[:] = basalFriction[:]
      netCDFfile.close()

      if not options.format_only:

#       Run CISM (NOTE two options here.)
        print options.parallel
        print 'Running',options.executable,'for experiment',experiment.upper(),'with domain size',size,'km'
        if options.parallel != None:
           if options.parallel > 0:
              exitCode = os.system('mpirun -np ' + str(options.parallel) + ' ./simple_glide '+filename+'.config')
           else:
              print 'Number of processors specified for parallel run is <=0.  Skipping the running of the model.'
              exitCode = 0
        else:
           exitCode = os.system('echo '+filename+'.config'+' | '+options.executable)

        if exitCode == 0:
           try:
    #         Extract the output data for comparison to the other models

    #         NOTE: The script now assumes that uvel_icegrid & vvel_icegrid are ALWAYS present.
    #         Those fields containthe ice velocity computed at the upper right corner of each grid cell.
    #         They appear to be on the x1,y1 grid in their metadata but are actually on the x0,y0 grid.
    #         The additional row/column include the first halo value past ewn/nsn.
    #         That value is valid at both 0.0 and 1.0 on the nondimensional coordinate system.
    #         Matrix manipulations for each test case below are done to create a larger matrix that goes from 0.0 to 1.0, inclusive.
    #         NOTE: The cases below are only writing [x,y,u,v] to the text file.  This is the minimum needed to compare to other models.
    #         In the future, the other additional fields specified in section 4 of http://homepages.ulb.ac.be/~fpattyn/ismip/ismiphom.pdf
    #         can be added.  wvel and the stresses are on the x1,y1 grid, so they would need to be interpolated to the x0,y0 grid
    #         since we are using that as the coordinate system in the text files.

    #         Open the netCDF file that was written by CISM
              netCDFfile = NetCDFFile(filename+'.out.nc','r')

              # Make x,y position arrays that can be used by all test cases.
              #   Want x/y positions to include the periodic edge at both the beginning and end
              xx = netCDFfile.variables['x0'][:]/(1000.0*float(size))
              xx = np.concatenate(([0.0],xx,[1.0]))
              yy = netCDFfile.variables['y0'][:]/(1000.0*float(size))
              yy = np.concatenate(([0.0],yy,[1.0]))
              if experiment in ('b','d'):
                 yy = yy[len(yy)/2]  # for the 2-d experiments, just use the middle y-index


              # Figure out u,v since all experiments needs at least one of them (avoids duplicate code in each case below
              us = netCDFfile.variables['uvel_icegrid'][0,0,:,:] * netCDFfile.variables['uvel_icegrid'].scale_factor  # top level of first time
              us = np.concatenate( (us[:,-1:], us), axis=1)  # copy the column at x=1.0 to x=0.0
              us = np.concatenate( (us[-1:,:], us), axis=0)  # copy the row at y=1.0 to y=0.0
              vs = netCDFfile.variables['vvel_icegrid'][0,0,:,:] * netCDFfile.variables['vvel_icegrid'].scale_factor  # top level of first time
              vs = np.concatenate( (vs[:,-1:], vs), axis=1)  # copy the column at x=1.0 to x=0.0
              vs = np.concatenate( (vs[-1:,:], vs), axis=0)  # copy the row at y=1.0 to y=0.0
              ub = netCDFfile.variables['uvel_icegrid'][0,-1,:,:] * netCDFfile.variables['uvel_icegrid'].scale_factor  # bottom level of first time
              ub = np.concatenate( (ub[:,-1:], ub), axis=1)  # copy the column at x=1.0 to x=0.0
              ub = np.concatenate( (ub[-1:,:], ub), axis=0)  # copy the row at y=1.0 to y=0.0
              vb = netCDFfile.variables['vvel_icegrid'][0,-1,:,:] * netCDFfile.variables['vvel_icegrid'].scale_factor  # bottom level of first time
              vb = np.concatenate( (vb[:,-1:], vb), axis=1)  # copy the column at x=1.0 to x=0.0
              vb = np.concatenate( (vb[-1:,:], vb), axis=0)  # copy the row at y=1.0 to y=0.0
              nan = ub*np.NaN  # create a dummy matrix for uncalculated values.

              # make arrays of the variables needed for each experiment
              # the icegrid velocities have the periodic edge in the last x-position.  We also want it in the first x-position.
              # After building the 2-d array as needed for each variable from the raw file data, then build a list called 'data'.
              if experiment == 'a':
                #  This is supposed to be: [('uvel',0),('vvel',0),('wvel',0),('tau_xz',-1),('tau_yz',-1),[deltap]]
                data = (us, vs, nan, nan, nan, nan)
              elif experiment == 'b':
                #  This is supposed to be: uvel(0), wvel(0), tau_xz(-1), deltaP
                data = (us, nan, nan, nan)
              elif experiment == 'c':
                #  This is supposed to be: [uvel',0),('vvel',0),('wvel',0),('uvel',-1),('vvel',-1),('tau_xz',-1),('tau_yz',-1), deltap]
                data = (us, vs, nan, ub, vb, nan, nan, nan)
              elif experiment == 'd':
                #  This is supposed to be:  [('uvel',0),('wvel',0),('uvel',-1),('tau_xz',-1), deltap]
                data = (us, nan, nan, nan, nan)
              elif experiment == 'f':
    #            variables = [('usurf',None),('uvel',0),('vvel',0),('wvel',0)]
                data = (nan, us, vs, nan)

    #         Write a "standard" ISMIP-HOM file (example file name: "cis1a020.txt") in the "output" subdirectory 
              ISMIP_HOMfilename = os.path.join('output',options.prefix+experiment+'%03d'%size+'.txt')
              ISMIP_HOMfile = open(ISMIP_HOMfilename,'w')
              for i, x in enumerate(xx):
                  for j, y in enumerate(yy):
                      if experiment in ('a','c','f'):  # include x and y positions
                        ISMIP_HOMfile.write('\t'.join(map(str,[x,y]+[v[j,i] for (v) in data]))+'\n')
                      else:  # only include x position
                        ISMIP_HOMfile.write('\t'.join(map(str,[x]+[v[j,i] for (v) in data]))+'\n')
              ISMIP_HOMfile.close()
              netCDFfile.close()
           except:
              print 'Error: The CISM output file for experiment '+experiment+' at size '+str(size)+' could NOT be read/post-processed successfully!'
        else:
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



