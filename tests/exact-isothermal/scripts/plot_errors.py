#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot absolute dome and max errors


import Numeric, Scientific.IO.NetCDF
import PyGMT, PyCF,sys

def parse_title(title):
    """Parse title string."""

    t = title.split(',')

    exp_name = t[0][-1]
    solver = t[1].strip()
    dx = t[2].strip()[:-2]
    dt = t[3].strip()[:-1]

    return (exp_name,solver,dx,dt)

if __name__ == '__main__':

    innames = sys.argv[1:-1]
    outname = sys.argv[-1]

    # start plotting
    plot = PyGMT.Canvas(outname,size='A4')
    plot.defaults['LABEL_FONT_SIZE']='12p'
    plot.defaults['ANOT_FONT_SIZE']='10p'
    
    key_y=2.5
    ysize = 7.
    dy = 1.
    
    bigarea = PyGMT.AreaXY(plot,size=[30,30])

    area_dome = PyGMT.AutoXY(bigarea,pos=[0.,key_y+dy],size=[15,ysize],logx=True)
    area_dome.xlabel = '@~D@~x [km]'
    area_dome.ylabel = 'absolute dome error [m]'
    area_dome.axis = 'WeSn'

    area_max = PyGMT.AutoXY(bigarea,pos=[0.,ysize+2*dy+key_y],size=[15,ysize],logx=True)
    area_max.xlabel = '@~D@~x [km]'
    area_max.ylabel = 'absolute maximum error [m]'
    area_max.axis = 'Wesn'
    
    s_keyarea = PyGMT.KeyArea(bigarea,pos=[0.,-0.8],size=[5,key_y])
    s_keyarea.num=[1,7]

    e_keyarea = PyGMT.KeyArea(bigarea,pos=[5.,-0.8],size=[10.,key_y])
    e_keyarea.num=[2,7]

    plotted = {}

    i = 0
    styles = {'non-lin':'','lin':'to','ADI':'ta'}
    colours = ['255/0/0','0/255/0','0/0/255','0/255/255','255/0/255','255/255/0','127/0/0','0/127/0','0/0/127','0/127/127','127/0/127','127/127/0']

    runs = {}

    # load data
    for f in innames:
        cffile = PyCF.CFloadfile(f)
        (exp_name,solver,dx,dt) = parse_title(cffile.title)
        if solver not in runs:
            runs[solver] = {}
        if exp_name not in runs[solver]:
            runs[solver][exp_name]={}
        diff = cffile.file.variables['thke'][-1,:,:] - cffile.file.variables['thk'][-1,:,:]
        centre = (Numeric.shape(diff)[0]-1)/2
        dome_e = diff[centre,centre]
        diff = Numeric.ravel(diff)
        max_e = max(abs(diff))
        runs[solver][exp_name][int(dx)] = [dome_e,max_e]
        cffile.close()

    # plot data
    done_grid = {}
    i = 0
    for s in runs:
        s_keyarea.plot_line(s,'3/0/0/0%s'%styles[s])
        for e in runs[s]:
            if e not in done_grid:
                done_grid[e] = i
                i = i + 1
                e_keyarea.plot_line('test %s'%e,'3/%s'%(colours[done_grid[e]]))

            dx = runs[s][e].keys()
            dx.sort()

            # get errors
            error_dome = []
            error_max = []
            for x in dx:
                error_dome.append(runs[s][e][x][0])
                error_max.append(runs[s][e][x][1])

            # plot line
            area_dome.line('-W3/%s%s'%(colours[done_grid[e]],styles[s]),dx,error_dome)
            area_dome.plotsymbol(dx,error_dome,size=0.1,args='-W1/%s'%(colours[done_grid[e]]))

            area_max.line('-W3/%s%s'%(colours[done_grid[e]],styles[s]),dx,error_max)
            area_max.plotsymbol(dx,error_max,size=0.1,args='-W1/%s'%(colours[done_grid[e]]))

    area_dome.finalise()
    area_dome.coordsystem()

    area_max.finalise()
    area_max.coordsystem()    

    plot.close()
