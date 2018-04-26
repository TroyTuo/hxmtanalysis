import numpy as np
import os
import pyfits as pf
from function import genlc
import matplotlib.pyplot as plt

det_box = [0,1,2]
pathlist = open('/home/tuoyl/Data_analysis/TimingTool/Glitch/'+\
        'le_pre_glitch_path.lst_test','r')
print pathlist
for path in pathlist:
    path = path[0:-1]
    print path
    if path == '':continue
    filename = path + 'le_screen_ang.fits'
    le_dir = path

    # eventype selection
    print "select eventype"
    hdulist = pf.open(filename)
    time = hdulist[1].data.field('time')
    et = hdulist[1].data.field('event_type')
    hdulist.close()
    det = det_box
    detnum = len(det)
    etindex = np.equal(et,1)
    ettime = time[etindex]
    etbin = np.arange(time[0],time[-1],1)
    et1rate,et1time = np.histogram(ettime,bins=etbin)
    et1time = et1time[:-1]
    et1time_tmp = et1time[:-1]
    evtindex = np.greater_equal(et1rate,detnum*1000)
    edgeindex = evtindex[1:]-evtindex[:-1]
    edges = et1time_tmp[edgeindex]
    if evtindex[0]:edges=np.concatenate([[time[0]],edges])
    if evtindex[-1]:edges=np.concatenate([edges,[time[-1]]])
    low_edge_tmp = edges[::2]
    up_edge_tmp = edges[1::2]
    low_edge = [x for x,y in zip(low_edge_tmp,up_edge_tmp) if y-x >= 2]
    up_edge = [y for x,y in zip(low_edge_tmp,up_edge_tmp) if y-x >= 2]

    x,y = genlc(time,pannel=False) 
    plt.subplot(2,1,1)
    plt.plot(x,y)
    x,y = genlc(time[etindex],pannel=False)
    plt.subplot(2,1,2)
    plt.plot(x,y)
    f = open(le_dir+'evttype.tim','w')
    for i in xrange(len(low_edge)):
        if i==len(low_edge)-1:
            wrt_str = '(time >= %s && time <= %s) '%(str(low_edge[i]),str(up_edge[i]))
        else:
            wrt_str = '(time >= %s && time <= %s) ||'%(str(low_edge[i]),str(up_edge[i]))
        f.write(wrt_str)
        print wrt_str
    f.close()
    
    ## select evttype=1 events
    fselect_text = 'fselect '+filename+' '+le_dir+'le_screen_evt1.fits @'+le_dir+'evttype.tim clobber=yes'
    print fselect_text
    os.system(fselect_text)
plt.show()
