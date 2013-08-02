#!/usr/bin/env python
# This script runs a "Stream / Shelf Experiment, based on the Goldberg et al. test case
# from JGR, 117, 2012, doi:10.1029/2011JF002246.

import sys, os, glob, shutil, numpy
from netCDF import *
from ConfigParser import ConfigParser
from math import sin, pi

# Check to see if a config file was specified on the command line.
# If not, confined-shelf.config is used.
if len(sys.argv) > 1:
  if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
    print '\nUsage:  python confined-shelf.py [FILE.CONFIG]\n'
    sys.exit(0)
  else:
    configfile = sys.argv[1]
else:
  configfile = 'stream-shelf.config'

# Check to see if #procs specified, relevant when running the code in parallel.
# If not, serial run (#procs==1) is performed. To run in parallel, the configure
# file must be specifed, but the nu,ber of processors does not
if len(sys.argv) > 2:
    nprocs = sys.argv[2]
else:
  nprocs = '1'

# Create a netCDF file according to the information in the config file.
parser = ConfigParser()
parser.read(configfile)
nx = int(parser.get('grid','ewn'))
ny = int(parser.get('grid','nsn'))
nz = int(parser.get('grid','upn'))
dx = float(parser.get('grid','dew'))
dy = float(parser.get('grid','dns'))
filename = parser.get('CF input', 'name')

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

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x.tolist()
netCDFfile.createVariable('y1','f',('y1',))[:] = y.tolist()

netCDFfile.createVariable('x0','f',('x0',))[:] = (dx/2 + x[:-1]).tolist()
netCDFfile.createVariable('y0','f',('y0',))[:] = (dy/2 + y[:-1]).tolist()

# *SFP* this has been changed so that the default value for 'flwa' is the same 
# as in the EISMINT-shelf test documentation, tests 3 & 4, found at:
# http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html

# Calculate values for the required variables.
thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg  = numpy.zeros([1,ny,nx],dtype='float32')
beta = numpy.empty([1,ny-1,nx-1],dtype='float32')
kbc  = numpy.zeros([1,ny-1,nx-1],dtype='int')
acab = numpy.zeros([1,ny,nx],dtype='float32') # *sfp* added acab field for prog. runs 
zero = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

# init vel fields (need to be zeroed at margins to set up kinbc of u=v=0 in seacism)
uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

# feel out other fields
acab[:] = 0.3           # value from Goldberg et al. (2012) test case (Table 1 - m/yr)
beta[0,:,:] = 30.4582   # value from Goldberg et al. (2012) test case (Table 1 - converted from Pa s/m to Pa yr/m)

# Domain size
Lx = (nx-1)*dx
Ly = (ny-1)*dy

# define bed topog
for i in range(nx):
  xx = float(i) * dx #- 0.5   # -1/2 < x < 1/2 
  for j in range(ny):
    yy = float(j) * dy #- 0.5 # -1/2 < y < 1/2
#    topg[0,j,i] = -(300.0 + 600.0 * sin( pi*(xx-2.0*dx) / (Lx-4.0*dx) ) )    # shelf front at domain bottom
    topg[0,j,i] = -(300.0 + 600.0 * sin( pi*(yy-2.0*dy) / (Ly-4.0*dy) ) )    # shelf front at domain right 

# shelf front at domain bottom 
#thk[0,4:-2,2:-2] = 500.     
#kbc[0,ny-4:,:]  = 1
#kbc[0,:,:3] = 1
#kbc[0,:,nx-4:] = 1
#topg[0,:,:4] = -300.0
#topg[0,:,nx-4:] = -300.0
#acab[0,ny-3:,:]  = 0    # zero out accum at edges to avoid buildup where u=0
#acab[0,:,:3] = 0
#acab[0,:,nx-3:] = 0

# shelf front at domain right     
thk[0,2:-2,2:-4] = 500.     
kbc[0,:,:3]  = 1
kbc[0,:3,:] = 1
kbc[0,ny-4:,:] = 1
topg[0,:4,:] = -300.0
topg[0,ny-4:,:] = -300.0
acab[0,:,:3]  = 0    # zero out accum at edges to avoid buildup where u=0
acab[0,:3,:] = 0
acab[0,ny-3:,:] = 0

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk',      'f',('time','y1','x1'))[:] = thk.tolist()
netCDFfile.createVariable('acab',     'f',('time','y1','x1'))[:] = acab.tolist()
netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kbc.tolist()
netCDFfile.createVariable('topg',     'f',('time','y1','x1'))[:] = topg.tolist()
netCDFfile.createVariable('beta',     'f',('time','y0','x0'))[:] = beta.tolist()
netCDFfile.createVariable('uvel',  'f',('time','level','y0','x0'))[:] = zero.tolist()

# *sfp* first option below adds ice stream vel profile for kin bc at upstream end
# *sfp* ... comment out for standard test case
#netCDFfile.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = vvel.tolist()
netCDFfile.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = zero.tolist()

netCDFfile.close()

# Run Glimmer
print 'Running Glimmer/CISM'
if len(sys.argv) > 2:
#   os.system('aprun -n'+nprocs+' ./simple_glide '+configfile+'')  # support for MPI runs on Jaguar
   os.system('mpirun -np '+nprocs+' ./simple_glide '+configfile+'')  # support for MPI runs on other machines
else:
   os.system('echo '+configfile+' | simple_glide')

# Clean up by moving extra files written by Glimmer to the "scratch" subdirectory
# Look for files with extension "txt", "log", or "nc"
for files in glob.glob('*.txt')+glob.glob('*.log'):
# Delete any files already in scratch with these filenames 
  if files in os.listdir('scratch'):
    os.remove(os.path.join('scratch',files))
# Move the new files to scratch
  shutil.move(files,'scratch')
