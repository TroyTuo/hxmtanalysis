#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
A script for generating phase-coherent timing models (ephemerides)
based on the rotating frequency searches and fittings with the raw
ephemerides as inputs.
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits
import sys, os
from scipy.optimize import curve_fit

def ccf(f1,f2):
    '''f1 is the original signal
       f2 is probe signal(shift and test)'''
    f2range = xrange(len(f2))
    mean_f1 = np.mean(f1)
    mean_f2 = np.mean(f2)
    delta_f1 = f1 - mean_f1
    delta_f2 = f2 - mean_f2
    sigma_f1 = np.sqrt(np.sum([x**2 for x in f1]))
    sigma_f2 = np.sqrt(np.sum([x**2 for x in f2]))
    y = [ np.sum(delta_f1 * np.roll(delta_f2,x))/(sigma_f1 * sigma_f2) for x in xrange(len(f2)) ]
    delay = np.where(y==max(y))[0]
    return y,delay

#read datafile
def read_file(datafile, colname='TDB', time_cut=30):
    '''time_cut is ephemeris effective time range(in unit of days)'''
    print "...cutting data into %g-day intervals..."%time_cut
    time = []
    hdulist = fits.open(datafile)
    time_all = hdulist[1].data.field(colname) #TODO:cut every 30 days?
    if time_cut == 0 or (max(time_all)-min(time_all))/86400.0 < time_cut:
        time.append(time_all)
    else:
        edges = np.arange(min(time_all),max(time_all),time_cut*86400)
        edges = np.array([5.19180059e+08, 5.21772059e+08, 5.24364059e+08, 5.26956059e+08,
            5.29548059e+08, 5.32140059e+08, 5.34732059e+08, 5.37324059e+08,
            5.39916059e+08, 5.43888000e+08, 5.45100059e+08, 5.47692059e+08,
            5.50284059e+08, 5.52876059e+08, 5.55468059e+08])
        for i in xrange(len(edges[0:-1])):
            edge0 = edges[i]
            edge1 = edges[i+1]
            data = time_all[(time_all>=edge0) & (time_all<=edge1)]
            time.append(data)
    return time
    

#read par file
def read_par(parname):
    pardata = open(parname,'r')
    stdpar = []
    parameters = np.zeros(13)
    for par in pardata:
        par = par[0:(len(par)-1)]
        stdpar.append(par)
    pardata.close()
    for i in xrange(len(stdpar)):
        if stdpar[i][:6]=='PEPOCH':
            PEPOCH_lst = stdpar[i].split(' ');PEPOCH = [x for x in PEPOCH_lst if x is not ''][1]
            parameters[0] = np.longdouble(PEPOCH) 
        if stdpar[i][:2]=='F0': 
            F0_lst = stdpar[i].split(' ');F0 = [x for x in F0_lst if x is not ''][1]
            parameters[1] = np.longdouble(F0) 
        if stdpar[i][:2]=='F1':
            F1_lst = stdpar[i].split(' ');F1 = [x for x in F1_lst if x is not ''][1]
            parameters[2] = np.longdouble(F1)
        if stdpar[i][:2]=='F2':
            F2_lst = stdpar[i].split(' ');F2 = [x for x in F2_lst if x is not ''][1]
            parameters[3] = np.longdouble(F2)
        if stdpar[i][:2]=='F3':
            F3_lst = stdpar[i].split(' ');F3 = [x for x in F3_lst if x is not ''][1]
            parameters[4] = np.longdouble(F3)
        if stdpar[i][:2]=='F4':
            F4_lst = stdpar[i].split(' ');F4 = [x for x in F4_lst if x is not ''][1]
            parameters[5] = np.longdouble(F4)
        if stdpar[i][:2]=='F5':
            F5_lst = stdpar[i].split(' ');F5 = [x for x in F5_lst if x is not ''][1]
            parameters[6] = np.longdouble(F5)
        if stdpar[i][:2]=='F6':
            F6_lst = stdpar[i].split(' ');F6 = [x for x in F6_lst if x is not ''][1]
            parameters[7] = np.longdouble(F6)
        if stdpar[i][:2]=='F7':
            F7_lst = stdpar[i].split(' ');F7 = [x for x in F7_lst if x is not ''][1]
            parameters[8] = np.longdouble(F7)
        if stdpar[i][:2]=='F8':
            F8_lst = stdpar[i].split(' ');F8 = [x for x in F8_lst if x is not ''][1]
            parameters[9] = np.longdouble(F8)
        if stdpar[i][:2]=='F9':
            F9_lst = stdpar[i].split(' ');F9 = [x for x in F9_lst if x is not ''][1]
            parameters[10] = np.longdouble(F9)
        if stdpar[i][:5]=='START':
            START_lst = stdpar[i].split(' ');START = [x for x in START_lst if x is not ''][1]
            parameters[11] = np.longdouble(START) 
        if stdpar[i][:6]=='FINISH':
            FINISH_lst = stdpar[i].split(' ');FINISH = [x for x in FINISH_lst if x is not ''][1]
            parameters[12] = np.longdouble(FINISH) 

    print "...finish reading ephemeris file..."
    return parameters

def write_par(inparname, outparname, EPOCH=0, f0=0, f1=0, f2=0, tstart=0, tfinish=0):
    with open(inparname, 'r') as f:
        lines = f.readlines()
    with open(outparname, 'w') as fout:
        for i in xrange(len(lines)):
            if lines[i][:6] == 'PEPOCH':
                lines[i] = 'PEPOCH      '+str(EPOCH)+' 1 0 \n'
            if lines[i][:5] == 'START':
                lines[i] = 'START      '+str(tstart)+' \n'
            if lines[i][:6] == 'FINISH':
                lines[i] = 'FINISH      '+str(tfinish)+'\n'
            if lines[i][:2] == 'F0':
                lines[i] = 'F0          '+str(f0)+' 1 0 \n'
            if lines[i][:2] == 'F1':
                lines[i] = 'F1          '+str(f1)+' 1 0 \n'
            if lines[i][:2] == 'F2':
                lines[i] = 'F2          '+str(f2)+' 1 0 \n'
            if lines[i][:2] == 'P0':
                #lines[i] = 'P0          '+str(1/f0)+' 1 0 \n'
                lines[i] = '\n'
            if lines[i][:2] == 'P1':
                #lines[i] = 'P1          '+str(-f1/f0**2)+' 1 0 \n'
                lines[i] = '\n'
            #TODO: modify P, P1, P2 as well?
            fout.write(lines[i])
    return
    

def fsearch(time, parfile, duration, fstep=0, frange=0,fig_flag=False,
        bin_profile=1000, mission='hxmt', threshold=0.16):
    '''
    Search frequency based on Pearson Chi-square test.
    The reference time is the start time of each time intervals.
    '''

    #set reference MJD 
    if mission=='hxmt' or mission=='HXMT':
        print '...mission is HXMT...'
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif mission=='fermi' or mission=='FERMI':
        print '...mission is FERMI...'
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910
    else:
        raise IOError('mission "%s" not supported, add refenece MJD data manually'%mission)

    #read parfile and parameters
    parameters = read_par(parfile)
    PEPOCH = parameters[0]
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    F0 = parameters[1]
    F1 = parameters[2]
    F2 = parameters[3]
    F3 = parameters[4]
    F4 = parameters[5]
    F5 = parameters[6]
    F6 = parameters[7]
    F7 = parameters[8]
    F8 = parameters[9]
    F9 = parameters[10]
    T_START = parameters[11]
    T_FINISH = parameters[12]

    fbest_all = np.array([])
    tf_all = np.array([])
    ferr_all = np.array([])

    # seperate data
    if duration == 0 or duration >= (max(time)-min(time)):
        edges = np.array([min(time),max(time)])
    else:
        edges = np.arange(min(time),max(time),duration)

    print "...searching frequency..."
    for i in xrange(len(edges[0:-1])):
        edge0 = edges[i]
        edge1 = edges[i+1]
        
        data = time[ (time>=edge0) & (time<=edge1) ]
        if len(data)==0:continue
        t0 = min(data) # search frequency, reference time is edge of time interval
        T0 = t0/86400.0 + MJDREFF + MJDREFI
        dt = t0 - pepoch 
        f0 = F0 + F1*dt + (1/np.math.factorial(2))*F2*(dt**2) + (1/np.math.factorial(3))*F3*(dt**3) + (1/np.math.factorial(4))*F4*(dt**4) + (1/np.math.factorial(5)*F5*(dt**5)) + (1/np.math.factorial(6)*F6*(dt**6)) + (1/np.math.factorial(7)*F7*(dt**7)) + (1/np.math.factorial(8)*F8*(dt**8)) + (1/np.math.factorial(9)*F9*(dt**9)) 
        f1 = F1 + F2*dt + (1/np.math.factorial(2))*F3*(dt**2) + (1/np.math.factorial(3))*F4*(dt**3) + (1/np.math.factorial(4))*F5*(dt**4) + (1/np.math.factorial(5)*F6*(dt**5)) + (1/np.math.factorial(6)*F7*(dt**6)) + (1/np.math.factorial(7)*F8*(dt**7)) + (1/np.math.factorial(8)*F9*(dt**8))
        f2 = F2 + F3*dt + (1/np.math.factorial(2))*F4*(dt**2) + (1/np.math.factorial(3))*F5*(dt**3) + (1/np.math.factorial(4))*F6*(dt**4) + (1/np.math.factorial(5)*F7*(dt**5)) + (1/np.math.factorial(6)*F8*(dt**6)) + (1/np.math.factorial(7)*F9*(dt**7))
        f3 = F3 + F4*dt + (1/np.math.factorial(2))*F5*(dt**2) + (1/np.math.factorial(3))*F6*(dt**3) + (1/np.math.factorial(4))*F7*(dt**4) + (1/np.math.factorial(5)*F8*(dt**5)) + (1/np.math.factorial(6)*F9*(dt**6))
        f4 = F4 + F5*dt + (1/np.math.factorial(2))*F6*(dt**2) + (1/np.math.factorial(3))*F7*(dt**3) + (1/np.math.factorial(4))*F8*(dt**4) + (1/np.math.factorial(5)*F9*(dt**5))
        f5 = F5 + F6*dt + (1/np.math.factorial(2))*F7*(dt**2) + (1/np.math.factorial(3))*F8*(dt**3) + (1/np.math.factorial(4))*F9*(dt**4)
        f6 = F6 + F7*dt + (1/np.math.factorial(2))*F8*(dt**2) + (1/np.math.factorial(3))*F9*(dt**3)
        f7 = F7 + F8*dt + (1/np.math.factorial(2))*F9*(dt**2)
        f8 = F8 + F9*dt
        f9 = F9 

        f = np.arange(f0-frange,f0+frange,fstep)
        chi_square = np.array([])
        N = len(data)
        b = N/bin_profile
        #for j in tqdm(range(0,len(f))):
        for j in xrange(len(f)):
            #phi_tmp = np.mod((data-t0)*f[j] ,1.0)
            phi_tmp = np.mod((data-t0)*f[j] + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4 +
                    (1/np.math.factorial(6))*((data-t0)**6)*f5 + (1/np.math.factorial(7))*((data-t0)**7)*f6 + (1/np.math.factorial(8))*((data-t0)**8)*f7 + 
                    (1/np.math.factorial(9))*((data-t0)**9)*f8 + (1/np.math.factorial(10))*((data-t0)**10)*f9 ,1.0)
            p_num = np.histogram(phi_tmp,bin_profile)[0]
            chi_square = np.append(chi_square,(np.std(p_num)**2/np.mean(p_num)))

        if max(chi_square)==min(chi_square):print('skip this profile');continue


        fbest = f[np.where(chi_square==max(chi_square))][0]
        fbest_all = np.append(fbest_all, fbest)
        tf_all = np.append(tf_all, T0)
        ferr_all = np.append(ferr_all, fstep)

        if fig_flag:
            plt.figure()
            plt.plot(f,chi_square)

    return fbest_all, tf_all, ferr_all

def toa_cal(time, parfile, duration=0, bin_profile=50, std_pro_file='', fig_flag=False, mission='HXMT',threshold=0.16):
    if mission=='hxmt' or mission=='HXMT':
        print '...mission is HXMT...'
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif mission=='fermi' or mission=='FERMI':
        print '...mission is FERMI...'
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910
    else:
        raise IOError('mission "%s" not supported, add refenece MJD data manually'%mission)

    #read parfile and parameters
    parameters = read_par(parfile)
    PEPOCH = parameters[0]
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    F0 = parameters[1]
    F1 = parameters[2]
    F2 = parameters[3]
    F3 = parameters[4]
    F4 = parameters[5]
    F5 = parameters[6]
    F6 = parameters[7]
    F7 = parameters[8]
    F8 = parameters[9]
    F9 = parameters[10]
    T_START = parameters[11]
    T_FINISH = parameters[12]

    # seperate data
    if duration == 0 or duration >= (max(time)-min(time)):
        edges = np.array([min(time),max(time)])
    else:
        edges = np.arange(min(time),max(time),duration)
    print "...calculating ToAs, every "+str(duration)+" seconds..."

    toa = np.array([])
    # use two method
    for i in xrange(len(edges[0:-1])):
        edge0 = edges[i]
        edge1 = edges[i+1]
        data = time[ (time>=edge0) & (time<=edge1) ]
        if len(data)==0:continue
        t0 = min(data)
        #t0 = pepoch #fold profile, reference time is PEPOCH
        T0 = t0/86400 + MJDREFF + MJDREFI
        dt = t0 - pepoch 
        f0 = F0 + F1*dt + (1/np.math.factorial(2))*F2*(dt**2) + (1/np.math.factorial(3))*F3*(dt**3) + (1/np.math.factorial(4))*F4*(dt**4) + (1/np.math.factorial(5)*F5*(dt**5)) + (1/np.math.factorial(6)*F6*(dt**6)) + (1/np.math.factorial(7)*F7*(dt**7)) + (1/np.math.factorial(8)*F8*(dt**8)) + (1/np.math.factorial(9)*F9*(dt**9)) 
        f1 = F1 + F2*dt + (1/np.math.factorial(2))*F3*(dt**2) + (1/np.math.factorial(3))*F4*(dt**3) + (1/np.math.factorial(4))*F5*(dt**4) + (1/np.math.factorial(5)*F6*(dt**5)) + (1/np.math.factorial(6)*F7*(dt**6)) + (1/np.math.factorial(7)*F8*(dt**7)) + (1/np.math.factorial(8)*F9*(dt**8))
        f2 = F2 + F3*dt + (1/np.math.factorial(2))*F4*(dt**2) + (1/np.math.factorial(3))*F5*(dt**3) + (1/np.math.factorial(4))*F6*(dt**4) + (1/np.math.factorial(5)*F7*(dt**5)) + (1/np.math.factorial(6)*F8*(dt**6)) + (1/np.math.factorial(7)*F9*(dt**7))
        f3 = F3 + F4*dt + (1/np.math.factorial(2))*F5*(dt**2) + (1/np.math.factorial(3))*F6*(dt**3) + (1/np.math.factorial(4))*F7*(dt**4) + (1/np.math.factorial(5)*F8*(dt**5)) + (1/np.math.factorial(6)*F9*(dt**6))
        f4 = F4 + F5*dt + (1/np.math.factorial(2))*F6*(dt**2) + (1/np.math.factorial(3))*F7*(dt**3) + (1/np.math.factorial(4))*F8*(dt**4) + (1/np.math.factorial(5)*F9*(dt**5))
        f5 = F5 + F6*dt + (1/np.math.factorial(2))*F7*(dt**2) + (1/np.math.factorial(3))*F8*(dt**3) + (1/np.math.factorial(4))*F9*(dt**4)
        f6 = F6 + F7*dt + (1/np.math.factorial(2))*F8*(dt**2) + (1/np.math.factorial(3))*F9*(dt**3)
        f7 = F7 + F8*dt + (1/np.math.factorial(2))*F9*(dt**2)
        f8 = F8 + F9*dt
        f9 = F9 
    
        phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4 +
                (1/np.math.factorial(6))*((data-t0)**6)*f5 + (1/np.math.factorial(7))*((data-t0)**7)*f6 + (1/np.math.factorial(8))*((data-t0)**8)*f7 + 
                (1/np.math.factorial(9))*((data-t0)**9)*f8 + (1/np.math.factorial(10))*((data-t0)**10)*f9 ,1.0)
        bin_x = np.arange(0,1,1.0/bin_profile)
        bin_x = np.append(bin_x,1.0)
        p_num_unnorm = np.histogram(phi,bin_x)[0]
        #p_num_unnorm = np.histogram(phi,bin_profile)[0]
        p_num = p_num_unnorm
        p_num = [(x - min(p_num))/(max(p_num)-min(p_num)) for x in p_num] # Normalization
        p_num_x = np.arange(0,1,1.0/bin_profile)
        
        for x in p_num:
            print x
        
        if np.std(p_num,ddof=1)<threshold:print('skip this profile');continue
        if std_pro_file != '':
            # standard profile
            std_pro = np.loadtxt(std_pro_file)
            std_pro = [(x - min(std_pro))/(max(std_pro)-min(std_pro)) for x in std_pro] # Normalization
        else:
            try:std_pro
            except:std_pro = p_num
        ## ccf shift
        y, delay = ccf(p_num,std_pro)
        p_num_std = np.roll(std_pro,delay)
        phi_peak = p_num_x[np.where(p_num_std==max(p_num_std))][0]

        if fig_flag:
            plt.figure()
            plt.step(p_num_x,p_num)
            plt.plot(p_num_x,p_num_std,color='red')

        toa_t0 = T0
        toa_f0 = f0
        toa_tmp = toa_t0 + (1/toa_f0) * phi_peak/86400
        toa = np.append(toa, toa_tmp)
    return toa

def write_toa(filename,toa,instru,error=0):
    with open(filename,'w')as f:
        f.write('FORMAT 1\n')
        for i in xrange(len(toa)):
            wrt_str = instru+'.toa 3000 %.20f %g bat\n'%(toa[i],error)
            f.write(wrt_str)


def ffitting(mjd, frequency):
    def fun(dT, f0, f1, f2):
        y = f0 + f1*dT + (1.0/2.0)*f2*dT**2
        return y
    t0 = mjd[0]
    dT = mjd-t0
    popt, pcov = curve_fit(fun, dT, frequency)
    print "...fitting frequency..."
    f0 = popt[0]
    f1 = popt[1]/86400.0
    f2 = popt[2]/(86400.0**2)
    return t0, f0, f1, f2

def save_frequency(outfrefile, time, frequency, ferror):
    with open(outfrefile, 'w')as fout:
        for i in xrange(len(time)):
            wrt_str = '%f %f %f \n'%(time[i], frequency[i], ferror[i])
            fout.write(wrt_str)
    return

def pass_argument():
    pass

def main(time, parfile, outparfile, outtimfile, outfrefile, time_cut, threshold=0.16):
    #fsearch
    fbest, mjd, ferr = fsearch(time, parfile, duration=3600*10, fig_flag=False,
            fstep=1e-7, frange=1e-5, bin_profile=50, mission='fermi',threshold=threshold)
    #save frequency 
    save_frequency(outfrefile, mjd, fbest, ferr)
    EPOCH, f0, f1, f2 = ffitting(mjd, fbest)
    #generate new parfile
    tstart = np.floor(min(mjd))
    tfinish = np.ceil(max(mjd+time_cut))
    #NOTE:after the first round of making parfile, comment write_par, update ToAs and parfile NOTE#
    #write_par(parfile, outparfile, EPOCH=EPOCH, f0=f0, f1=f1, f2=f2, tstart=tstart, tfinish=tfinish)
    #calculate ToAs
    toa = toa_cal(time, outparfile, duration=3600*20, bin_profile=50, std_pro_file='test_profile.dat', fig_flag=False, mission='FERMI',threshold=threshold)
    write_toa(outtimfile, toa, 'Fermi', error=2000)
    os.system('/Users/tuoyouli/Documents/tempo2/bin/tempo2 -gr plk -f %s %s'%(outparfile, outtimfile))
    #raw_input("before continuing, refit the Ephemeris using TEMPO2, Press the <ENTER> key to continue...")

    return outparfile


if __name__ == '__main__':
    datafile = '../geminga_filtered_gti_bary.fits' #TODO:KeyboardInput
    #rawparfile = '../atnf_test.par'
    rawparfile = '../results/geminga_eph_0.par'
    newparfile_prefix = '../results/geminga_eph'
    newtimfile_prefix = '../results/geminga_toa'
    newfrefile_prefix = '../results/geminga_fre'
    mission = 'FERMI'
    time_cut = 30
    time_list = read_file(datafile, colname='Time', time_cut=time_cut) 

    tmp_i = 0
    for time, i in zip(time_list, xrange(len(time_list)))[tmp_i:]:
        if len(time) == 0:continue
        print "...processing %g/%g time interval..."%(i, len(time_list)-1)
        outparfile = newparfile_prefix+'_%g.par'%i
        outtimfile = newtimfile_prefix+'_%g.tim'%i
        outfrefile = newfrefile_prefix+'_%g.dat'%i
        if i == tmp_i:
            lastparfile = main(time, rawparfile, outparfile, outtimfile, outfrefile, time_cut)
        else:
            lastparfile = main(time, lastparfile,  outparfile, outtimfile, outfrefile, time_cut)
        plt.show()

