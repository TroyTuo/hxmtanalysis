#!/usr/bin/env python
'''Dead Time correction for HXMT payloads'''
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
from function import genlc

def exp_cal(time):
    #calculate exposure time
    tmp_time = np.roll(time,-1);tmp_time[-1]=time[-1]
    time_step = tmp_time - time
    unexp_flag = np.greater(time_step,1)
    unexp_dur = np.sum(time_step[unexp_flag])
    exposure = max(time)-min(time)-unexp_dur
    return exposure

def deadtime_cal(time_det_num,time_dtime, DeadTime_PHODet_num):
    ''' calculate deadtime and exposure with 
    time_det_num(arrival time for specific detector),
    time_dtime(time of deadtime file), 
    DeadTime_PHODet(deadtime data in deadtime file) as inputs.
    And return deadtime, exposure time'''
    ######## find the edge time of light curve for on detector ########
    time = time_det_num; DeadTime_PHODet = DeadTime_PHODet_num;
    #time = np.sort(time)
    print 'no sort'
    # upper limits
    tmp_time = np.roll(time,-1);tmp_time[-1]=time[-1]
    time_step = tmp_time - time
    edge_flag = np.greater(time_step,1)
    edges = time[edge_flag]
    # lowwer limits
    tmp_time = np.roll(time,1);tmp_time[0]=time[0]
    time_step = tmp_time - time
    edge_flag = np.less(time_step,-1)
    edges = np.append(edges, time[edge_flag])
    # min & max limits
    edges = np.append(edges, min(time))
    edges = np.append(edges, max(time))
    edges = np.sort(edges)
    left_edges = edges[::2]
    right_edges = edges[1::2]
#temporary
    time_dtime_sel=np.array([])
    time_dtime_sel1=np.array([])
    DeadTime_PHODet_sel = np.array([])
    DeadTime_PHODet_sel1 = np.array([])
    plt.subplot(2,1,2)
    
    for i in xrange(len(left_edges)):
        index_0 = int(np.where(time_dtime==np.ceil(left_edges[i]))[0])
        index_1 = int(np.where(time_dtime==np.ceil(right_edges[i]))[0])
        time_dtime_sel = np.append(time_dtime_sel, time_dtime[index_0:index_1])
        time_dtime_sel1 = np.append(time_dtime_sel1, np.mean(time_dtime[index_0:index_1]))
        DeadTime_PHODet_sel = np.append(DeadTime_PHODet_sel, DeadTime_PHODet[index_0:index_1])
        DeadTime_PHODet_sel1 = np.append(DeadTime_PHODet_sel1, np.mean(DeadTime_PHODet[index_0:index_1]))
    plt.plot(time_dtime_sel,DeadTime_PHODet_sel)
    plt.plot(time_dtime_sel1,DeadTime_PHODet_sel1,'*-',markersize=8)
#temporary
    
    ######## find the time intervals in dead time file #########
    index = np.array([]); deadtime = 0
    for i in edges:
        index = np.append(index, np.where(time_dtime==np.floor(i))[0][0])
    lower_index = index[::2]
    upper_index = index[1::2]
    # dead time in lower edge intervals
    for i in xrange(len(lower_index)):
        lower_deadtime = (np.floor(left_edges[i]) - left_edges[i]) * DeadTime_PHODet[int(lower_index[i])]
    # dead time in upper edge intervals
    for i in xrange(len(upper_index)):
        upper_deadtime = (right_edges[i] - np.ceil(right_edges[i])) * DeadTime_PHODet[int(upper_index[i])]
    # dead time of other interval 
    for i in xrange(len(lower_index)):
        deadtime = deadtime + np.sum(DeadTime_PHODet[int(lower_index[i])+1:int(upper_index[i])-1])
        deadtime = deadtime + upper_deadtime - lower_deadtime
    exposure = exp_cal(time)
    print "deadtime, exposure, counts rate",deadtime, exposure, len(time)/exposure
    return deadtime, exposure


def read_deadtimefile(filename):
    ######### read deadtime files ########
    dtimefile = filename
    filelist = open(dtimefile,'r')
    time_dtime = np.array([]); 
    DeadTime_PHODet_0 = np.array([])
    DeadTime_PHODet_1 = np.array([])
    DeadTime_PHODet_2 = np.array([])
    DeadTime_PHODet_3 = np.array([])
    DeadTime_PHODet_4 = np.array([])
    DeadTime_PHODet_5 = np.array([])
    DeadTime_PHODet_6 = np.array([])
    DeadTime_PHODet_7 = np.array([])
    DeadTime_PHODet_8 = np.array([])
    DeadTime_PHODet_9 = np.array([])
    DeadTime_PHODet_10 = np.array([])
    DeadTime_PHODet_11 = np.array([])
    DeadTime_PHODet_12 = np.array([])
    DeadTime_PHODet_13 = np.array([])
    DeadTime_PHODet_14 = np.array([])
    DeadTime_PHODet_15 = np.array([])
    DeadTime_PHODet_16 = np.array([])
    DeadTime_PHODet_17 = np.array([])
    for file in filelist:
        file = file[0:-1]
        if file == '':continue
        print file
        hdulist = pf.open(file)
        tb = hdulist[1].data
        time_dtime = np.append(time_dtime,tb.field('Time'))
        DeadTime_PHODet_0 = np.append(DeadTime_PHODet_0,tb.field('DeadTime_PHODet_0')/1e6)
        DeadTime_PHODet_1 = np.append(DeadTime_PHODet_1,tb.field('DeadTime_PHODet_1')/1e6)
        DeadTime_PHODet_2 = np.append(DeadTime_PHODet_2,tb.field('DeadTime_PHODet_2')/1e6)
        DeadTime_PHODet_3 = np.append(DeadTime_PHODet_3,tb.field('DeadTime_PHODet_3')/1e6)
        DeadTime_PHODet_4 = np.append(DeadTime_PHODet_4,tb.field('DeadTime_PHODet_4')/1e6)
        DeadTime_PHODet_5 = np.append(DeadTime_PHODet_5,tb.field('DeadTime_PHODet_5')/1e6)
        DeadTime_PHODet_6 = np.append(DeadTime_PHODet_6,tb.field('DeadTime_PHODet_6')/1e6)
        DeadTime_PHODet_7 = np.append(DeadTime_PHODet_7,tb.field('DeadTime_PHODet_7')/1e6)
        DeadTime_PHODet_8 = np.append(DeadTime_PHODet_8,tb.field('DeadTime_PHODet_8')/1e6)
        DeadTime_PHODet_9 = np.append(DeadTime_PHODet_9,tb.field('DeadTime_PHODet_9')/1e6)
        DeadTime_PHODet_10 = np.append(DeadTime_PHODet_10,tb.field('DeadTime_PHODet_10')/1e6)
        DeadTime_PHODet_11 = np.append(DeadTime_PHODet_11,tb.field('DeadTime_PHODet_11')/1e6)
        DeadTime_PHODet_12 = np.append(DeadTime_PHODet_12,tb.field('DeadTime_PHODet_12')/1e6)
        DeadTime_PHODet_13 = np.append(DeadTime_PHODet_13,tb.field('DeadTime_PHODet_13')/1e6)
        DeadTime_PHODet_14 = np.append(DeadTime_PHODet_14,tb.field('DeadTime_PHODet_14')/1e6)
        DeadTime_PHODet_15 = np.append(DeadTime_PHODet_15,tb.field('DeadTime_PHODet_15')/1e6)
        DeadTime_PHODet_16 = np.append(DeadTime_PHODet_16,tb.field('DeadTime_PHODet_16')/1e6)
        DeadTime_PHODet_17 = np.append(DeadTime_PHODet_17,tb.field('DeadTime_PHODet_17')/1e6)
    print "Finish reading deadtime files"
    return time_dtime, DeadTime_PHODet_0 ,DeadTime_PHODet_1,DeadTime_PHODet_2,DeadTime_PHODet_3,\
            DeadTime_PHODet_4,DeadTime_PHODet_5,DeadTime_PHODet_6,DeadTime_PHODet_7,DeadTime_PHODet_8,\
            DeadTime_PHODet_9,DeadTime_PHODet_10,DeadTime_PHODet_11,DeadTime_PHODet_12,DeadTime_PHODet_13,\
            DeadTime_PHODet_14,DeadTime_PHODet_15,DeadTime_PHODet_16,DeadTime_PHODet_17
