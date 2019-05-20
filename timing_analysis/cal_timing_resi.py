#!/usr/bin/env python

from __future__ import division
import numpy as np
from toa_module.fun import read_par
from toa_module.fun import phi_cal
import matplotlib.pyplot as plt


def read_tim(timfile):
    '''read .tim file'''

    toas = np.array([])

    filedata = open(timfile)
    lines = filedata.readlines()
    for i in xrange(len(lines)):
        if i == 0:continue
        line = lines[i]
        line = line[0:-1].split()
        if line == []:continue
        toa = line[2]
        error = float(line[3])
        toas = np.append(toas, toa)
    return toas, error*1e-6

def timing_resi(parfile,timfile):
    # read parfile
    print 'parfile: ',parfile
    print 'toafile: ',timfile
    parameters = read_par(parfile)
    toas, error0 = read_tim(timfile)
    toas = [float(x) for x in toas]

    PEPOCH = parameters[0]
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
    dt = (toas - PEPOCH)*86400

    # calculate phi variation
    phi = np.mod(((toas-PEPOCH)*86400)*F0 + (1/2)*(((toas-PEPOCH)*86400)**2)*F1 + (1/6)*(((toas-PEPOCH)*86400)**3)*F2 + (1/24)*(((toas-PEPOCH)*86400)**4)*F3 + (1/120)*(((toas-PEPOCH)*86400)**5)*F4 +
            (1/np.math.factorial(6))*(((toas-PEPOCH)*86400)**6)*F5 + (1/np.math.factorial(7))*(((toas-PEPOCH)*86400)**7)*F6 + (1/np.math.factorial(8))*(((toas-PEPOCH)*86400)**8)*F7 + 
            (1/np.math.factorial(9))*(((toas-PEPOCH)*86400)**9)*F8 + (1/np.math.factorial(10))*(((toas-PEPOCH)*86400)**10)*F9 ,1.0)
    residual = phi - np.mean(phi)

    # cal error in period
    f_real = F0 + F1*dt + (1/np.math.factorial(2))*F2*(dt**2) + (1/np.math.factorial(3))*F3*(dt**3) + (1/np.math.factorial(4))*F4*(dt**4) + (1/np.math.factorial(5)*F5*(dt**5)) + (1/np.math.factorial(6)*F6*(dt**6)) + (1/np.math.factorial(7)*F7*(dt**7)) + (1/np.math.factorial(8)*F8*(dt**8)) + (1/np.math.factorial(9)*F9*(dt**9)) 
    P_real = 1/f_real

    error = error0/P_real

    plt.plot(toas, residual,'.', marker='^',markersize=8,color='b')
    for i in xrange(len(toas)):
        plt.plot([toas[i],toas[i]],[residual[i]-error[i], residual[i]+error[i]],'-',color='b')
    plt.show()
    return toas, residual, error



#parfile = 'hxmt_pre_glitch.par'
#timfile = 'hxmt_pre_glitch.tim'
#toas0, residual0, error0 = timing_resi(parfile,timfile)
#with open('./results/resi_pre.dat','w')as f:
#    for i in xrange(len(toas0)):
#        wrt_str = '%f %f %f\n'%(toas0[i],residual0[i],error0[i])
#        f.write(wrt_str)
#
#parfile = 'hxmt_post_glitch.par'
#timfile = 'hxmt_post_glitch.tim'
#toas1, residual1, error1 = timing_resi(parfile,timfile)
#with open('results/resi_post1.dat','w')as f:
#    for i in xrange(len(toas1)):
#        wrt_str = '%f %f %f\n'%(toas1[i],residual1[i],error1[i])
#        f.write(wrt_str)

parfile = '/Users/tuoyouli/Work/Crab_analysis/Glitch/hxmt_ephemeris/newpar_1.par'
timfile = '/Users/tuoyouli/Work/Crab_analysis/Glitch/hxmt_ephemeris/hxmt_post_glitch.tim'
toas1, residual1, error1 = timing_resi(parfile,timfile)
with open('results/resi_post_new1.dat','w')as f:
    for i in xrange(len(toas1)):
        wrt_str = '%f %f %f\n'%(toas1[i],residual1[i],error1[i])
        f.write(wrt_str)
