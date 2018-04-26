from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
from fun import *
from peaksearch_ccf import peak_search
from scipy.optimize import curve_fit
from deadtime_cor import read_deadtimefile
from projection_cor import projection_cor

def profile_fitting(pathlistname, filename, parfilename, std_pro, dtimefile, f3_flag=True, f4_flag=True, fig_flag=False,dtime_flag=True,psf_flag=True):
    '''return epoch, rate and error of rate'''
    ########### generate clean profile with pulsed-photons only ##############
    profile_x = np.arange(0,0.999,1000)
    profile_x = np.linspace(0,0.999,1000)
    binsize = 1000
    counts = np.array([]); rate = np.array([]); epoch = np.array([]); error = np.array([])
    dtime_frac = np.array([])
    
    #  select photons in bk phases
    index_max = peak_search(std_pro, std_pro)[0]
    index_max = index_max[0]
    if index_max >= binsize:
        index_max = index_max-binsize
    phase_min = profile_x[index_max]+0.6
    phase_max = profile_x[index_max]+0.8 
    print "min_phase = ",phase_min,"max_phase=",phase_max
    
    #select those phase intervals
    phase_min = np.mod(phase_min,1)
    phase_max = np.mod(phase_max,1)
    if phase_min <= phase_max:
        std_bk = std_pro[np.logical_and(np.greater_equal(profile_x,phase_min),np.less_equal(profile_x,phase_max))] 
        std_bk_x = profile_x[np.logical_and(np.greater_equal(profile_x,phase_min),np.less_equal(profile_x,phase_max))]
    else:
        std_bk = std_pro[np.logical_or(np.less_equal(profile_x,phase_max),np.greater_equal(profile_x,phase_min))]
        std_bk_x = profile_x[np.logical_or(np.less_equal(profile_x,phase_max),np.greater_equal(profile_x,phase_min))]
    
    # exclude background photons for std profile
    std_pro = std_pro - np.mean(std_bk)
    std_pro = std_pro/np.sum(std_pro)
    counts_std_pro = np.sum(std_pro)
    plt.figure()
    plt.title('std profile')
    line1,=plt.plot(profile_x,std_pro,color='r')
    plt.axvline(x=phase_max,linewidth=0.5,color='black')
    plt.axvline(x=phase_min,linewidth=0.5,color='black')
    plt.axhline(y=0,linewidth=0.5,color='grey')
    ##############
    
    ############## Fitting profile by strandard profile ################
    
    ## define fitting function
    def func(x, a, b):
        return a * x + b
    tdb_all = np.array([]); time_all = np.array([]) 
    pathlist = open(pathlistname,'r')
    for path in pathlist:
        path = path[0:-1]
        if path == '':continue
        hdulist = pf.open(path + filename)
        tb = hdulist[1].data
        try:
            tdb = tb.field('TDB')
            time = tb.field('Time')
            detid = tb.field('DetID')
            tdb_all = np.append(tdb_all,tdb)
            time_all = np.append(time_all,time)
            print path
        except:
            print "NO TDB column"
            continue 
        hdulist.close()

    ######## read dead files ########
    time_dtime_tmp,\
    DeadTime_PHODet_0 ,DeadTime_PHODet_1,DeadTime_PHODet_2,DeadTime_PHODet_3,DeadTime_PHODet_4,DeadTime_PHODet_5,\
    DeadTime_PHODet_6,DeadTime_PHODet_7,DeadTime_PHODet_8,DeadTime_PHODet_9,DeadTime_PHODet_10,DeadTime_PHODet_11,\
    DeadTime_PHODet_12,DeadTime_PHODet_13,DeadTime_PHODet_14,DeadTime_PHODet_15,DeadTime_PHODet_16,DeadTime_PHODet_17 = read_deadtimefile(dtimefile)
    
    DeadTime_PHODet_tmp = (DeadTime_PHODet_0 +DeadTime_PHODet_1+DeadTime_PHODet_2+DeadTime_PHODet_3+DeadTime_PHODet_4+DeadTime_PHODet_5+\
    DeadTime_PHODet_6+DeadTime_PHODet_7+DeadTime_PHODet_8+DeadTime_PHODet_9+DeadTime_PHODet_10+DeadTime_PHODet_11+\
    DeadTime_PHODet_12+DeadTime_PHODet_13+DeadTime_PHODet_14+DeadTime_PHODet_15+DeadTime_PHODet_16+DeadTime_PHODet_17)/18

    time_range = np.array([min(time_all),179539200,180835200,184636800])
    time_range = np.arange(min(time_all),max(time_all),10000)

    ######## read psf file #########
    time_qt_all, psf3_all = projection_cor(pathlistname)

    for time_range_index in xrange(len(time_range)-1):
        #select time range
        tdb = tdb_all[ np.logical_and(time_all>=time_range[time_range_index], time_all<=time_range[time_range_index+1]) ]
        time = time_all[ np.logical_and(time_all>=time_range[time_range_index], time_all<=time_range[time_range_index+1]) ]
        if len(tdb)<=100:continue
    
        # fold profile for time rage
        profile_x,profile,toa_tmp = pfold(tdb,parfilename,duration=3600*99999,f3_flag=f3_flag,f4_flag=f4_flag,fig_flag=False,bin_profile=binsize,std_pro_file=std_pro,threshold=0)
        # fitting each profile by standar profile
        y, delay = ccf(profile, std_pro)
        std_pro = np.roll(std_pro, delay)
        popt, pcov = curve_fit(func, std_pro, profile)
        a = popt[0]; b = popt[1]
        perr = np.sqrt(np.diag(pcov))
        std_a = perr[0]; std_b = perr[1]
        #calculate exposure time
        exposure = exp_cal(time)

        #find edges
        left_edges, right_edges = find_edges(time)

        time_dtime = np.array([]);DeadTime_PHODet = np.array([])

        cor_factor = 1
        if dtime_flag:
            # select good time intervals
            for i in xrange(len(left_edges)):
                time_dtime = np.append(time_dtime, time_dtime_tmp[np.logical_and(time_dtime_tmp>=left_edges[i], time_dtime_tmp<=right_edges[i])])
                DeadTime_PHODet = np.append(DeadTime_PHODet, DeadTime_PHODet_tmp[np.logical_and(time_dtime_tmp>=left_edges[i], time_dtime_tmp<=right_edges[i])])
            dtime_cor = np.mean(DeadTime_PHODet)
            cor_factor = cor_factor*(1-1.8*dtime_cor)
        
        time_qt = np.array([]); psf = np.array([])
        if psf_flag:
            # select good time intervals
            for i in xrange(len(left_edges)):
                time_qt = np.append(time_qt, time_qt_all[np.logical_and(time_qt_all>=left_edges[i], time_qt_all<=right_edges[i])])
                psf = np.append(psf, psf3_all[np.logical_and(time_qt_all>=left_edges[i], time_qt_all<=right_edges[i])])
            psf_cor = np.mean(psf)
            cor_factor = cor_factor*psf_cor
        ######### Calculate counts rate of pulse photons ###############
        counts = np.append(counts, popt[0] * counts_std_pro/cor_factor)
        rate = np.append(rate, popt[0]*counts_std_pro/(exposure*cor_factor))
        error = np.append(error, std_a/a * rate)
        epoch = np.append(epoch, (min(time)+max(time))/2)




    

        if fig_flag:
            plt.figure()
            plt.plot(profile_x,profile,color='r')
            plt.plot(profile_x,popt[0]*std_pro+popt[1])
            # annotation text
            text = 'a*x+b \n'+'a= '+str(popt[0])+'\n b = '+str(popt[1])+'\n std_a = '+str(perr[0]) +'\n std_b ='+str(perr[1])
            plt.annotate(text, xy=(1, 1), xytext=(-15, -15), fontsize=10,
            xycoords='axes fraction', textcoords='offset points',
            bbox=dict(facecolor='white', alpha=0.8),
            horizontalalignment='right', verticalalignment='top')
    ##############
    
#    ############## Calculate counts rate of pulse photons ###############
#        counts = np.append(counts, popt[0] * counts_std_pro)
#        rate = np.append(rate, popt[0]*counts_std_pro/exposure)
#        error = np.append(error, std_a/a * rate)
#        epoch = np.append(epoch, (min(time)+max(time))/2) 

    print epoch, rate, error
    return epoch, rate, error


def find_edges(time):
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
    return left_edges, right_edges

def exp_cal(time):
    #calculate exposure time
    tmp_time = np.roll(time,-1);tmp_time[-1]=time[-1]
    time_step = tmp_time - time
    unexp_flag = np.greater(time_step,1)
    unexp_dur = np.sum(time_step[unexp_flag])
    exposure = max(time)-min(time)-unexp_dur
    return exposure









