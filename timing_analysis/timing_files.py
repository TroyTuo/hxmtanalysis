#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
A script for calculating frequency and folding profiles based on
phase-coherent timing models(ephemerisds)
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import fits
import glob
from natsort import natsorted
from function import ccf

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
def read_file(datafile, colname='TDB', parfile='', tstart=0, tstop=0, instrument='hxmt'):
    '''read data file with TSTART and TFINISH in parfile'''
    if instrument == 'hxmt' or instrument == 'HXMT':
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif instrument == 'fermi' or instrument == 'FERMI':
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910

    hdulist = fits.open(datafile)
    time = hdulist[1].data.field(colname) 
    if parfile:
        parameters = read_par(parfile)
        tstart = (parameters[11]-MJDREFF-MJDREFI)*86400 
        tstop  = (parameters[12]-MJDREFF-MJDREFI)*86400
        time = time[(time>=tstart) & (time<=tstop)]
    elif tstart != 0 or tstop != 0:
        time = time[time>=tstart & time<=tstop]
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

#read ToA file
def read_toa(timname):
    #with some bizarre lines
    file = open(timname)
    lines = file.readlines()
    file.close()
    toa = np.array([])
    for line in lines[1:]:
        line = line.split(' ')
        line = [x for x in line if x!='']
        if line[0] == 'C' :continue
        toa = np.append(toa, line[2])
    toa = np.asarray(toa, dtype=np.float128)

    #data = np.loadtxt(timname, skiprows=1, dtype={
    #    'names':('instru','what1','toa','err','what2'),
    #    'formats':('S','g','float128','float128','S')})
    #toa = np.asarray([x[2] for x in data])
    return toa


def profile_cal(time, parfile, timfile, bin_profile=50, instrument='hxmt'):
    if instrument == 'hxmt' or instrument == 'HXMT':
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif instrument == 'fermi' or instrument == 'FERMI':
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910

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

    data = time
    t0 = pepoch # !!! set one reference point
    T0 = t0/86400 + MJDREFF + MJDREFI
    f0 = F0
    f1 = F1
    f2 = F2
    f3 = F3
    f4 = F4
    f5 = F5
    f6 = F6
    f7 = F7
    f8 = F8
    f9 = F9 
    #print 'periodic parameters',f0,f1,f2,f3,f4,f5,f6,f7,f8,f9

    #phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2,1.0)
    #phi = np.mod(np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2,1.0)+0.5, 1) -0.5

    #calculate mean phi
    toas = read_toa(timfile)
    phi0, phi_toa = mean_phi(toas, 57900, T0, F0=f0, F1=f1, F2=f2) 

    #phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 - phi0,1.0)
    phi = np.mod(np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 - phi0,1.0)+0.5, 1) -0.5

    profile = np.histogram(phi,bin_profile)[0]

    return phi, profile, phi0

def pass_argument():
    pass

def mean_phi(toas, epoch_new, epoch_old, F0, F1,F2=0,F3=0,F4=0,F5=0,F6=0,F7=0,F8=0,F9=0):
    '''modify the parameters to asign epoch_new
    and get the delta Phi of a profile refered to epoch_new'''
    epoch_new = epoch_old
    dt = (epoch_new - epoch_old)*86400
    f0 = F0 + F1*dt + (1/2)*F2*(dt**2)
    f1 = F1 + F2*dt 
    f2 = F2

    dt_toa = (toas-epoch_new) * 86400
    phi = np.mod((dt_toa)*f0 + (1/2)*((dt_toa)**2)*f1 + (1/6)*((dt_toa)**3)*f2,1.0)
    #phi = np.mod(np.mod((dt_toa)*f0 + (1/2)*((dt_toa)**2)*f1 + (1/6)*((dt_toa)**3)*f2,1.0)+0.5, 1) -0.5
    mean_phi = np.mean(phi)

    ##NOTE:We calculate the mean value of Phase without those points deviating the mean value!NOTE
    #phi_tmp = np.asarray([x for x in phi if np.abs(x-mean_phi)<0.1])
    #phi = phi_tmp
    ####

    mean_phi = np.mean(phi)
    plt.figure()
    plt.plot(phi,'.')

    return mean_phi, phi

if __name__ == '__main__':
    datafile = '../geminga_filtered_gti_bary.fits'
    parfile_list = natsorted(glob.glob('../results/geminga_eph_*.par'), key=lambda y:y.lower())
    timfile_list = natsorted(glob.glob('../results/geminga_toa_*.tim'), key=lambda y:y.lower())
    profile_all = []
    phi_all = np.array([])
    time_all = np.array([])
    for parfile, timfile in zip(parfile_list, timfile_list)[:]:
        print parfile, timfile
        time = read_file(datafile, colname='Time', parfile=parfile, tstart=0, tstop=0, instrument='fermi')
        phi, profile, phi0 = profile_cal(time, parfile, timfile, bin_profile=100, instrument='fermi')
        print phi0, phi
        phi_all = np.append(phi_all, phi)
        time_all = np.append(time_all, time)
        shift_index = int(phi0*len(profile))
        profile_all.append(np.roll(profile, -shift_index))
print time_all, phi_all
H, x, y = np.histogram2d(time_all, phi_all, bins=(100, 100))
print x
plt.figure('3d')
plt.imshow(H,aspect='auto',origin='lower')
plt.figure()
for h in H:
    plt.plot(np.append(h, h))

plt.show()

#    plt.figure()
#    plt.imshow(profile_all,aspect='auto',origin='lower')
#    plt.show()










