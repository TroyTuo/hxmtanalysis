# -*- coding: <encoding name> -*-
from __future__ import division
import numpy as np
import pyfits as pf
import sys
import matplotlib.pyplot as plt
from function import genlc
from function import ccf

import os

def read_par(parname):
    pardata = open(parname,'r')
    stdpar = []
    parameters = np.zeros(9,dtype='longdouble')
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
    print 'parameters=',parameters
    return parameters

def read_data(datalistname):
    data = np.array([])
    datalist = open(datalistname,'r')
    for file in datalist:
        file = file[0:-1]
        print file
        try:
            hdulist = pf.open(file)
        except:
            print "no such file"
        try:
            tdb = hdulist[1].data.field('tdb')
        except:
            print "NO TDB column"
            continue
        data = np.append(data, tdb)
    data.sort()
    return data

def pfold(time,parfile,duration,f0_flag=True,f1_flag=True,f2_flag=True,f3_flag=False,f4_flag=False,f5_flag=True,f6_flag=True,f7_flag=True,fig_flag=True,bin_profile=1000,threshold=0,std_pro_file='',gen_std_pro=False,out_std_pro_file='',mission='hxmt'):
    ''' NOTICE: this pfold function is different with pfold in function base in hxmtanalysis,
    which fold the profile without parameter file'''
    if mission=='hxmt' or mission=='HXMT':
        print 'mission is HXMT...'
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif mission=='fermi' or mission=='FERMI':
        print 'mission is FERMI...'
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910
    else:
        print 'NO such mission'
    #read parfile and parameters
    parameters = read_par(parfile)
    PEPOCH = parameters[0]
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    if f0_flag:
        F0 = parameters[1]
    else:
        F0 = 0
    if f1_flag:
        F1 = parameters[2]
    else:
        F1 = 0
    if f2_flag:
        F2 = parameters[3]
    else:
        F2 = 0
    if f3_flag:
        F3 = parameters[4]
    else:
        F3 = 0
    if f4_flag:
        F4 = parameters[5]
    else:
        F4 = 0
    if f5_flag:
        F5 = parameters[6]
    else:
        F5 = 0
    if f6_flag:
        F6 = parameters[7]
    else:
        F6 = 0
    if f7_flag:
        F7 = parameters[8]
    else:
        F7 = 0

    # standard profile
    #std_pro = np.loadtxt(std_pro_file)
    #std_pro = [(x - min(std_pro))/(max(std_pro)-min(std_pro)) for x in std_pro] # Normalization

    # seperate data
    edges = np.arange(min(time),max(time),duration)
    if duration >= (max(time)-min(time)):
        edges = np.array([min(time),max(time)])
    p_x = []
    profile = []
    profile_std = []
    toa = []
    for i in xrange(len(edges[0:-1])):
        edge0 = edges[i]
        edge1 = edges[i+1]
        data = time[ (time>=edge0) & (time<=edge1) ]
        if len(data)==0:
            print "EMPTY"
            continue
        t0 = min(data)
        t0 = pepoch
        T0 = t0/86400 + MJDREFF + MJDREFI
        dt = t0 - pepoch 
        f0 = F0# + F1*dt + (1/2)*F2*(dt**2) + (1/6)*F3*(dt**3) + (1/24)*F4*(dt**4)
        f1 = F1# + F2*dt + (1/2)*F3*(dt**2) + (1/6)*F4*(dt**3)
        f2 = F2# + F3*dt + (1/2)*F4*(dt**2)
        f3 = F3# + F4*dt
        f4 = F4
        print 'f0,f1,f2,f3,f4',f0,f1,f2,f3,f4

        phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4,1.0)
        bin_x = np.arange(0,1,1.0/bin_profile)
        bin_x = np.append(bin_x,1.0)
        p_num_unnorm = np.histogram(phi,bin_x)[0]
        #p_num_unnorm = np.histogram(phi,bin_profile)[0]
        p_num = p_num_unnorm
        p_num = [(x - min(p_num))/(max(p_num)-min(p_num)) for x in p_num] # Normalization
        p_num_x = np.arange(0,1,1.0/bin_profile)
        if np.std(p_num,ddof=1)<threshold:continue

        if not gen_std_pro:
            #toa calculation
            ## ccf shift
            #y, delay = ccf(p_num,std_pro)
            #p_num_std = np.roll(std_pro,delay)
            #phi_peak = p_num_x[np.where(p_num_std==max(p_num_std))][0]
            phi_peak = p_num_x[np.where(p_num == max(p_num))][0]
            toa_tmp = T0 + (1/f0) * phi_peak
            toa_tmp = T0 + (1/f0) * phi_peak/86400
            p_x = np.append(p_x,p_num_x)
            profile = np.append(profile,p_num_unnorm)
            toa.append(toa_tmp)

    if gen_std_pro:
        with open(out_std_pro_file,'w')as f:
            for pi in xrange(len(p_num_unnorm)):
                wrt_str = '%f \n'%(p_num_unnorm[pi])
                f.write(wrt_str)

    if fig_flag:
        plt.figure('profiles')
        plt.plot(p_num_x,p_num,'r')
        ## plot light curve
        plt.figure('light curve')
        x,y = genlc(time,binsize=1)
        plt.plot(x,y)
        for i in xrange(len(edges[0:-1])):
             plt.axvline(x=edges[i],linewidth=0.5,color='r')
             plt.axvline(x=edges[i+1],linewidth=0.5,color='b')

    return p_x,profile,toa

def fsearch(time,parfile,duration,fstep,frange,f0_flag=True,f1_flag=True,f2_flag=True,f3_flag=False,f4_flag=False,fig_flag=False,bin_cs=20,bin_profile=1000,threshold=3e6,mission='hxmt'):
    ''' NOTICE: this fsearch function is different with fsearch in function base in hxmtanalysis,
    which calculate best frequency without parameter file'''
    if mission=='hxmt' or mission=='HXMT':
        print 'mission is HXMT...'
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif mission=='fermi' or mission=='FERMI':
        print 'mission is FERMI...'
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910
    else:
        print 'NO such mission'
    #read parfile and parameters
    parameters = read_par(parfile)
    print 'after read',type(parameters[0]),type(parameters[1])
    PEPOCH = parameters[0]
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    print 'pepoch=',pepoch
    if f0_flag:
        F0 = parameters[1]
    else:
        F0 = 0
    if f1_flag:
        F1 = parameters[2]
    else:
        F1 = 0
    if f2_flag:
        F2 = parameters[3]
    else:
        F2 = 0
    if f3_flag:
        F3 = parameters[4]
    else:
        F3 = 0
    if f4_flag:
        F4 = parameters[5]
    else:
        F4 = 0

    edges = np.arange(min(time),max(time),duration)
    if duration >= (max(time)-min(time)):
        edges = np.array([min(time),max(time)])
    tf = np.zeros(len(edges[0:-1]),dtype='longdouble')
    fre = np.zeros(len(edges[0:-1]),dtype='longdouble')
    fre_err = np.zeros(len(edges[0:-1]),dtype='longdouble')
    for i in xrange(len(edges[0:-1])):
        edge0 = edges[i]
        edge1 = edges[i+1]
        
        data = time[ (time>=edge0) & (time<=edge1) ]
        if len(data)==0:continue
        if len(data)<=threshold:continue
        print min(data),max(data),len(data)
        t0 = min(data)
        T0 = t0/86400.0 + MJDREFF + MJDREFI
        dt = t0 - pepoch 
        print 'cal f0',type(F0)
        f0 = F0 + F1*dt + (1.0/2)*F2*(dt**2) + (1.0/6)*F3*(dt**3) + (1.0/24)*F4*(dt**4)
        f1 = F1 + F2*dt + (1.0/2)*F3*(dt**2) + (1.0/6)*F4*(dt**4)
        f2 = F2 + F3*dt + (11.0/2)*F4*(dt**2)
        f3 = F3 + F4*dt
        f4 = F4
        print 'f0,f1,f2=',
        print '{0:.15f}'.format(f0),
        print '{0:.15f}'.format(f1),f2,type(F0),type(f0)

        f = np.arange(f0-frange,f0+frange,fstep)
        print type(f[0])
        chi_square = np.zeros(len(f))
        chi_square = np.array([])
        N = len(data)
        b = N/bin_cs
        for f1 in np.arange(-3.72e-10,-3.69e-10,0.5e-12):
            for j in range(0,len(f)):
                phi_tmp = np.mod((data-t0)*f[j] + (1.0/2)*((data-t0)**2)*f1 + (1.0/6.0)*((data-t0)**3)*f2 + (1.0/24)*(data-t0)**f3,1.0)
                p_num = np.histogram(phi_tmp,bin_cs)[0]
                #chi_square[j] = np.std(p_num)**2/np.mean(p_num)
                chi_square = np.append(chi_square,(np.std(p_num)**2/np.mean(p_num)))
                with open('./f0-f1search/f0-f1_search'+str(j)+'.dat','w')as ffile:
                    for ichi in chi_square:
                        ffile.write('%f \n'%ichi)
            

            percent = float(j)*100/len(f)
            sys.stdout.write(" fsearch complete: %.2f"%percent);
            sys.stdout.write("%\r");
            sys.stdout.flush()
            print "f: ",f
            print "f1: ",np.arange(-3.72e-10,-3.69e-10,0.5e-12)
            #fbest = f[np.where(chi_square==max(chi_square))][0]
            #fre[i] = fbest
            #tf[i] = T0
            #fre_err[i] = fstep/(max(chi_square)/bin_cs)
            #if fig_flag:
            plt.figure()
            #plt.plot(f,chi_square)
            plt.plot(chi_square)
    tf = tf[tf>0]
    fre =fre[fre>0]
    fre_err = fre_err[fre_err>0]

#    return p_num_x_2,p_num_2,f,chi_square
    return tf,fre,fre_err

def toa_write(filename,toa,instru):
    with open(filename,'w')as f:
        f.write('FORMAT 1\n')
        for i in xrange(len(toa)):
            wrt_str = instru+'.toa 3000 %.20f 33 bat\n'%(toa[i])
            f.write(wrt_str)

def fre_write(filename,mjd, fre, fre_err):
    with open(filename,'w')as f:
        for i in xrange(len(fre)):
            wrt_str = '%.16f %.16f %.16f \n'%(mjd[i],fre[i],fre_err[i])
            f.write(wrt_str)

def lucy_next_res(prob, now_res, raw_res):

    count = len(now_res)
    temp = np.inner(prob, now_res)
    temp2 = np.inner(raw_res / temp, np.transpose(prob))
    next_res = now_res * temp2
    return next_res

def lucy_iteration(prob, raw_res, n):
    for i in range(n):
        if i == 0:
            res = lucy_next_res(prob, raw_res, raw_res)
        else:
            res = lucy_next_res(prob, res, raw_res)
    return res

def le_delay(tdb):
    delay_time_tmp = []
    while len(delay_time_tmp) < len(tdb):
        print 1
        x = np.random.uniform(0.2e-3,1.18e-3,2*len(tdb))
        y = np.random.uniform(0,0.02,2*len(tdb))
        print 2
        delay_time_tmp = x[y-20.408*x+4.081e-3 <= 0]
        print 3
    delay_time = random.sample(delay_time_tmp,len(tdb))
    delay_time = np.asarray(delay_time)
    print delay_time
    return delay_time

def nor_pro_err(y):
    '''calculation error of normalized profile'''
    y = np.asarray(y)
    sig_y = np.sqrt(y)
    sig_max = np.sqrt(max(y))
    sig_min = np.sqrt(min(y))
    sigma_square = (sig_y**2 + sig_min**2)/((max(y)-min(y))**2) + (((y - min(y))**2)*(sig_max**2 + sig_min**2))/((max(y) - min(y))**4)
    sigma = np.sqrt(sigma_square)
    return sigma


def phi_cal(time,parfile,f0_flag=True,f1_flag=True,f2_flag=True,f3_flag=False,f4_flag=False):
    MJDREFF = 0.0007660185
    MJDREFI = 55927
    #read parfile and parameters
    parameters = read_par(parfile)
    PEPOCH = parameters[0]
    pepoch = (PEPOCH - MJDREFF - MJDREFI)*86400
    if f0_flag:
        F0 = parameters[1]
    else:
        F0 = 0
    if f1_flag:
        F1 = parameters[2]
    else:
        F1 = 0
    if f2_flag:
        F2 = parameters[3]
    else:
        F2 = 0
    if f3_flag:
        F3 = parameters[4]
    else:
        F3 = 0
    if f4_flag:
        F4 = parameters[5]
    else:
        F4 = 0

    data = time
#    t0 = min(data)
    t0 = pepoch # !!! set one reference point
    T0 = t0/86400 + MJDREFF + MJDREFI
    dt = t0 - pepoch 
    f0 = F0# + F1*dt + (1/2)*F2*(dt**2) + (1/6)*F3*(dt**3) + (1/24)*F4*(dt**4)
    f1 = F1# + F2*dt + (1/2)*F3*(dt**2) + (1/6)*F4*(dt**3)
    f2 = F2# + F3*dt + (1/2)*F4*(dt**2)
    f3 = F3# + F4*dt
    f4 = F4
    print "periodic parameters: ",f0,f1,f2,f3,f4

    ## write column to fits file
    #calc_text = 'ftcalc '+infile+' '+outfile+' Phase '+\
    #        '"((TDB-'+str(t0)+')*'+'('+str(f0)+')'+\
    #        '+ (1/2)*((TDB-'+str(t0)+')**2)*'+'('+str(f1)+')'+\
    #        '+ (1/6)*((TDB-'+str(t0)+')**3)*'+'('+str(f2)+')'+\
    #        '+ (1/24)*((TDB-'+str(t0)+')**4)*'+'('+str(f3)+')'+\
    #        '+ (1/120)*((TDB-'+str(t0)+')**5)*'+'('+str(f4)+')'+\
    #        ')%1"'+' clobber=yes'
    #print calc_text
    #os.system(calc_text)
    phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4,1.0)
    return phi
