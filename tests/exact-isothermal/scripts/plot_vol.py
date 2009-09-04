#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot relative volume errors


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

def calc_verror(cffile):

    verror = (cffile.file.variables['ivol'][:]-cffile.file.variables['ivole'][:])/cffile.file.variables['ivole'][:]
    verror[0] = 0.
    return verror

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
    vol_err = PyGMT.AutoXY(bigarea,pos=[0.,key_y+dy],size=[15,ysize])
    vol_err.xlabel = 'time [ka]'
    vol_err.ylabel = 'relative volume error'
    vol_err.axis = 'WeSn'
    
    s_keyarea = PyGMT.KeyArea(bigarea,pos=[0.,-0.8],size=[5,key_y])
    s_keyarea.num=[1,7]

    e_keyarea = PyGMT.KeyArea(bigarea,pos=[5.,-0.8],size=[10.,key_y])
    e_keyarea.num=[2,7]

    plotted = {}

    i = 0
    styles = {'non-lin':'','lin':'to','ADI':'ta'}
    colours = ['255/0/0','0/255/0','0/0/255','0/255/255','255/0/255','255/255/0','127/0/0','0/127/0','0/0/127','0/127/127','127/0/127','127/127/0']
    for f in innames:
        cffile = PyCF.CFloadfile(f)
        (exp_name,solver,dx,dt) = parse_title(cffile.title)
        title = '@~D@~x=%skm, @~D@~t=%sa'%(dx,dt)
        if title not in plotted:
            plotted[title] = i
            e_keyarea.plot_line(title,'3/%s'%(colours[i]))
            i = i+1
        vol_err.line('-W3/%s%s'%(colours[plotted[title]],styles[solver]),cffile.time(None),calc_verror(cffile))
        cffile.close()

    for s in styles:
        s_keyarea.plot_line(s,'3/0/0/0%s'%styles[s])

    vol_err.finalise()
    vol_err.coordsystem()

    bigarea.text([7.5,ysize+2*dy+key_y],"Experiment %s"%exp_name,textargs='20 0 0 CM')
    
    plot.close()
