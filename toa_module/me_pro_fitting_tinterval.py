from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
from fun import *
from peaksearch_ccf import peak_search
from scipy.optimize import curve_fit
from projection_cor import projection_cor_le
from function import exp_cal
from function import find_edges

def me_profile_fitting_ti(pathlist, filename, parfilename, std_pro,duration=300, f3_flag=True, f4_flag=True, fig_flag=False,dtime_flag=False,psf_flag=False):
    '''return epoch, rate and error of rate'''
    ########### generate clean profile with pulsed-photons only ##############
    profile_x = np.arange(0,0.999,1000)
    profile_x = np.linspace(0,0.999,1000)
    binsize = 1000
    counts = np.array([]); rate = np.array([]); epoch = np.array([]); error = np.array([])
    dtime_frac = np.array([])
    frac_index = np.array([])
    
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
    for path in pathlist:
        path = path[0:-1]
        if path == '':continue
        try:
            hdulist = pf.open(path + filename)
        except:
            print "datafile not exist"
        tb = hdulist[1].data
        try:
            tdb_all = tb.field('TDB')
            time_all = tb.field('Time')
            channel = tb.field('Channel')
            print path
        except:
            print "NO TDB column"
            continue 
        hdulist.close()
    
        # fold profile for each time interval
        edges = np.arange(min(time_all),max(time_all),duration)
        for i in xrange(len(edges[0:-1])):
            tdb = tdb_all[(time_all>=edges[i]) & (time_all<=edges[i+1])]
            if len(tdb)==0:
                print 'empty',
                continue
            if max(tdb)-min(tdb)<=200:
                print 'empty',
                continue
            time = time_all[(time_all>=edges[i]) & (time_all<=edges[i+1])]
            profile_x,profile,toa_tmp = pfold(tdb,parfilename,duration=3600*99999,f3_flag=f3_flag,f4_flag=f4_flag,fig_flag=False,bin_profile=binsize,std_pro_file=std_pro,threshold=0)
            # fitting each profile by standar profile
            y, delay = ccf(profile, std_pro)
            std_pro = np.roll(std_pro, delay)
            popt, pcov = curve_fit(func, std_pro, profile)
            a = popt[0]; b = popt[1]
            perr = np.sqrt(np.diag(pcov))
            std_a = perr[0]; std_b = perr[1]
            #calculate exposure time
            exposure = exp_cal(tdb)

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
        
        ############## Calculate counts rate of pulse photons ###############
            counts = np.append(counts, popt[0] * counts_std_pro)
            rate = np.append(rate, popt[0]*counts_std_pro/exposure)
            error = np.append(error, std_a/a * rate)
            epoch = np.append(epoch, (min(time)+max(time))/2) 

        ############# psf correction ################
            if psf_flag:
                time_qt, psf1, psf2, psf3 = projection_cor_le(pathlist)
                psf3_frac = np.mean(psf3[(time_qt >= edges[i]) & (time_qt <= edges[i+1])])
                frac_index = np.append(frac_index, psf3_frac)

            if dtime_flag:
        ############## dead time correction ############### 
                hdulist = pf.open(path + 'me_dtime.fits')
                tb = hdulist[1].data
                time_dtime = tb.field('Time')
                deadtime = (tb.field('DEADTIME0')+\
                        tb.field('DEADTIME1')+\
                        tb.field('DEADTIME2')+\
                        tb.field('DEADTIME3')+\
                        tb.field('DEADTIME4')+\
                        tb.field('DEADTIME5')+\
                        tb.field('DEADTIME6')+\
                        tb.field('DEADTIME7')+\
                        tb.field('DEADTIME8'))/9
                left_edges,right_edges = find_edges(time)
                print left_edges,right_edges
                deadtime_i = np.array([])
                for i in xrange(len(left_edges)):
                    deadtime_i = np.append(deadtime_i,\
                            deadtime[(time_dtime>=left_edges[i])\
                            &(time_dtime<=right_edges[i])])
                print np.sum(deadtime_i)
                dtime_frac = np.append(dtime_frac, np.sum(deadtime_i)/exposure)
    print "dtime_frac: ",dtime_frac
    if dtime_flag:
        return epoch, rate, error, dtime_frac
    if psf_flag:
        return epoch, rate, error, frac_index
    else:
        return epoch, rate, error












