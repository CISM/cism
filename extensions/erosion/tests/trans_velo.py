#!/usr/bin/env python

# Copyright 2005, Magnus Hagdorn

# create basal velocity field for testing transport code

import math, PyCF, Numeric

def create_velo(outname, grid_spacing):

    # all lengths are in kilometres
    sizex = [0.,1500.]  # area extent in x
    sizey = [0.,1500.]  # area extent in y
    centrex = 50.       # x centre of segment
    centrey = 50.       # y centre of segment
    radius = 1400.      # radius of circle
    start = 30.         # start angle
    sweep = 30.         # sweep of circle
    sedx=[50.,150.]      # extent of initial sediment cover in x
    sedy=[50.,150.]      # extent of initial sediment cover in y
    velo=66.66           # magnitude of velocity m/a
    sed_thck=5.         # initial sediment thickness

    numx = int((sizex[1]-sizex[0])/grid_spacing)+1
    numy = int((sizey[1]-sizey[0])/grid_spacing)+1

    start = math.radians(start)
    end = start+math.radians(sweep)

    # create file, etc.
    cffile = PyCF.CFcreatefile(outname)
    cffile.title = "Transport Test %f km grid"%grid_spacing
    cffile.institution = "University of Edinburgh"

    cffile.createDimension('x0',numx-1)
    cffile.createDimension('x1',numx)
    cffile.createDimension('y0',numy-1)
    cffile.createDimension('y1',numy)
    cffile.createDimension('level',1)
    cffile.createDimension('time',None)
    #creating variables
    varx=cffile.createVariable('x1')
    varx[:] = ((sizex[0]+grid_spacing*Numeric.arange(numx))*1000.).astype(Numeric.Float32)
    varx=cffile.createVariable('x0')
    varx[:] = ((sizex[0]+grid_spacing/2.+grid_spacing*Numeric.arange(numx-1))*1000.).astype(Numeric.Float32)
    
    vary=cffile.createVariable('y1')
    vary[:] = ((sizey[0]+grid_spacing*Numeric.arange(numy))*1000.).astype(Numeric.Float32)
    vary=cffile.createVariable('y0')
    vary[:] = ((sizey[0]+grid_spacing/2.+grid_spacing*Numeric.arange(numy-1))*1000.).astype(Numeric.Float32)
    
    varlevel=cffile.createVariable('level')
    varlevel[0] = 1

    vartime=cffile.createVariable('time')
    vartime[0] = 0

    # to metres
    centrex = centrex*1000.
    centrey = centrey*1000.
    radius = radius*1000.
        
    # create velos
    ubas = cffile.createVariable('ubas')
    vbas = cffile.createVariable('vbas')
    ubas[0,:,:] = 0.
    vbas[0,:,:] = 0.

    ubas[0,:,:] = velo
    vbas[0,:,:] = velo

##    for j in range(0,numy-1):
##        y = vary[j]-centrex
##        for i in range(0,numx-1):
##            x = varx[i]-centrex
##            dist = math.sqrt( x**2 + y**2)
##            a=math.atan2(y,x)
##            if dist <= radius and a >=start and a <=end:
##                ubas[0,j,i] = math.cos(a)*velo
##                vbas[0,j,i] = math.sin(a)*velo


    ubas[0,0,:] = 0.
    ubas[0,-1,:] = 0.
    ubas[0,:,0] = 0.
    ubas[0,:,-1] = 0.
    vbas[0,0,:] = 0.
    vbas[0,-1,:] = 0.
    vbas[0,:,0] = 0.
    vbas[0,:,-1] = 0.

    cffile.close()

if __name__ == '__main__':

    create_velo("velo.20km.nc",20.)
    create_velo("velo.10km.nc",10.)
    create_velo("velo.5km.nc",5.)
