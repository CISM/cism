#!/usr/bin/env python
#
# Magnus Hagdorn
#
# plot difference between numeric and exact solution

import os.path
import PyCF,PyGMT, sys, Numeric


if __name__ == '__main__':

    parser = PyCF.CFOptParser()
    parser.width = 12.5
    parser.time()
    parser.plot()
    opts = PyCF.CFOptions(parser,2)

    infile = opts.cffile()
    time = opts.times(infile)

    exact = infile.getvar('thke')
    diff = infile.getvar('thk').getGMTgrid(time)
    diff.data = Numeric.transpose(infile.file.variables['thk'][time,:,:] - infile.file.variables['thke'][time,:,:])

    plot = opts.plot()
    plot.defaults['LABEL_FONT_SIZE']='12p'
    plot.defaults['ANOT_FONT_SIZE']='10p'
    bigarea = PyGMT.AreaXY(plot,size=[30,30])

    area = PyCF.CFArea(bigarea,infile,pos=[0.,3.],size=opts.options.width)
    area.raw_image(infile,time,diff,'../data/error.cpt')
    area.contour(exact,[0.01],'-W2/0/0/0',time)
    area.coordsystem()
    PyGMT.colourkey(area,'../data/error.cpt',title='H@-num@--H@-exact@-',pos=[0.,-1.75],size=[opts.options.width,.75],args='-L')

    bigarea.text([opts.options.width/2.,opts.options.width+4.],"Experiment %s"%infile.title,textargs='20 0 0 CM')

    plot.close()
