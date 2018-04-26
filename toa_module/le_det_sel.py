#!/usr/bin/env python
''' select LE detector ID'''
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import os

pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
        'le_pre_glitch_path.lst_test'
infile = 'le_screen_evt1.fits'
outfile = 'le_screen_sml.fits'
with open(pathlistname,'r')as f:
    lines = f.read()
    pathlist = lines.split('\n')
    pathlist = pathlist[0:-1]

detlist = np.loadtxt('/home/tuoyl/Downloads/le_detid.dat')
print detlist[0][0]
for j in xrange(len(pathlist)):
    with open(pathlist[j] + 'det_sel.tim','w')as f:
        for i in xrange(len(detlist)):
            if detlist[i][1] == 1.0:
                print detlist[i][0]
                wrt_str = 'Det_ID == '+str(detlist[i][0]) + ' ||'
            f.write(wrt_str)
        print wrt_str
#for i in xrange(len(pathlist)):
#    print pathlist[i]

