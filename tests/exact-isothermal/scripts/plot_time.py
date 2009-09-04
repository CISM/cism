#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot time steps and run time


import os.path
import PyGMT, sys

def parse_title(title):
    """Parse title string."""

    config = title.split('_')
    exp = config[1]
    solver = config[2]
    dx = int(config[3][:-2])
    dt = float(config[4][:-1])

    return (exp,solver,dx,dt)

if __name__ == '__main__':

    inname = sys.argv[1]
    outname = sys.argv[2]

    # parse results file
    runs = {}
    for line in file(inname).readlines():
        # ignore comments and empty lines
        line = line.strip()
        pos = line.find('#')
        if pos>-1:
            line = line[:pos]
        if len(line)==0:
            continue
        line = line.split()
        fname = line[0][:-len('.config')]
        (exp_name,solver,dx,dt) = parse_title(fname)
        time = float(line[1])
        if solver not in runs:
            runs[solver] = {}
        if exp_name not in runs[solver]:
            runs[solver][exp_name]={}
        runs[solver][exp_name][dx] = [time,dt]

    # start plotting
    plot = PyGMT.Canvas(outname,size='A4')
    plot.defaults['LABEL_FONT_SIZE']='12p'
    plot.defaults['ANOT_FONT_SIZE']='10p'
    
    key_y=2.5
    ysize = 7.
    dy = 1.
    
    bigarea = PyGMT.AreaXY(plot,size=[30,30])

    area_dt = PyGMT.AutoXY(bigarea,pos=[0.,key_y+dy],size=[15,ysize],logx=True,logy=True)
    area_dt.xlabel = '@~D@~x [km]'
    area_dt.ylabel = '@~D@~t [a]'
    area_dt.axis = 'WeSn'

    area_rt = PyGMT.AutoXY(bigarea,pos=[0.,ysize+2*dy+key_y],size=[15,ysize],logx=True,logy=True)
    area_rt.xlabel = '@~D@~x [km]'
    area_rt.ylabel = 'run time [s]'
    area_rt.axis = 'Wesn'
    
    s_keyarea = PyGMT.KeyArea(bigarea,pos=[0.,-0.8],size=[5,key_y])
    s_keyarea.num=[1,7]

    e_keyarea = PyGMT.KeyArea(bigarea,pos=[5.,-0.8],size=[10.,key_y])
    e_keyarea.num=[2,7]

    styles = {'non-lin':'','lin':'to','ADI':'ta'}
    colours = ['255/0/0','0/255/0','0/0/255','0/255/255','255/0/255','255/255/0','127/0/0','0/127/0','0/0/127','0/127/127','127/0/127','127/127/0']

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

            rt = []
            dt = []
            for x in dx:
                rt.append(runs[s][e][x][0])
                dt.append(runs[s][e][x][1])

            # plot line
            area_rt.line('-W3/%s%s'%(colours[done_grid[e]],styles[s]),dx,rt)
            area_rt.plotsymbol(dx,rt,size=0.1,args='-W1/%s'%(colours[done_grid[e]]))

            area_dt.line('-W3/%s%s'%(colours[done_grid[e]],styles[s]),dx,dt)
            area_dt.plotsymbol(dx,dt,size=0.1,args='-W1/%s'%(colours[done_grid[e]]))            
    
    area_dt.finalise()
    area_dt.coordsystem()

    area_rt.finalise()
    area_rt.coordsystem()    

    plot.close()
