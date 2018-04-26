#!/usr/bin/env python
''' exclude slew data by RA, DEC of HXMT pointing direction'''
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
import os

pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/me_post_glitch_path.lst'
#pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/he_post_glitch_path.lst'
#pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/he_pre_glitch_path.lst_test'
#pathlistname = '/home/tuoyl/Data_analysis/TimingTool/Glitch/le_post_glitch_path.lst'
#filename = 'he_screen_NaI.fits'
filename = 'me_screen.fits'
#filename = 'he_screen_NaI.fits'
#filename = 'le_screen.fits'
#filename = 'le_screen.fits'
#outname = 'he_screen_NaI_ang.fits'
outname = 'me_screen_ang.fits'
#outname = 'he_screen_NaI_ang.fits'
#outname = 'le_screen_ang.fits'
#outname = 'le_screen_ang.fits'
#outname = 'le_screen_ang.fits'
ra_crab = 83.63322083
dec_crab = 22.01446111111 
ang_criterion = 0.1


with open(pathlistname,'r')as f:
    lines = f.read()
    pathlist = lines.split('\n')
    pathlist = pathlist[0:-1]

plt.figure()
fignum = np.ceil(np.sqrt(len(pathlist)))
print fignum
for i in xrange(len(pathlist)):
    ehkfile = pathlist[i]+'/../AUX/EHK.fits'
    try:
        hdulist = pf.open(ehkfile)
    except:
        print "######## EHK ERROR ##########",pathlist[i]
        continue
    tb = hdulist[1].data
    time_ehk = tb.field('time')
    hdulist.close()
    ra = tb.field('ra')
    dec = tb.field('dec')

    ang_dist = np.sqrt( (ra-ra_crab)**2 + (dec-dec_crab)**2 )
    print ra,dec
    print ang_dist
    index = np.less_equal(ang_dist,ang_criterion)
    index_tmp = np.append(False,index)
    index_tmp = np.append(index_tmp,False)
    edgeindex = index_tmp[1:]-index_tmp[:-1]
    time_ehk_tmp = np.append(time_ehk[0],time_ehk)
    edges = time_ehk_tmp[edgeindex]
    low_edges = edges[::2]
    up_edges = edges[1::2]
    with open(pathlist[i]+'/ang_sel.tim','w')as f:
        for j in xrange(len(low_edges)):
            if j == len(low_edges)-1:
                wrt_str =  '(time > %s && time < %s) '%(str(low_edges[j]),str(up_edges[j]))
            else:
                wrt_str = '(time > %s && time < %s) ||'%(str(low_edges[j]),str(up_edges[j]))
            f.write(wrt_str)
            print wrt_str
    plt.subplot(fignum,fignum,i+1)
    plt.plot(time_ehk,ra)
    plt.plot(time_ehk,dec)
    for j in xrange(len(low_edges)):
        plt.axvline(x=low_edges[j],color='r')
        plt.axvline(x=up_edges[j],color='b')

    fselect_text = 'fselect ' + pathlist[i] + filename + ' ' + pathlist[i] + outname + ' @'+pathlist[i]+'/ang_sel.tim clobber=yes > $glitch/ang_sel.log'
    print fselect_text
    os.system(fselect_text)
plt.show()





