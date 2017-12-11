# -*- coding: <encoding name> -*-
from __future__ import division
import numpy as np
import pyfits as pf
import sys
import matplotlib.pyplot as plt

def read_par(parname):
    pardata = open(parname,'r')
    stdpar = []
    parameters = np.zeros(6,dtype='longdouble')
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
            print F0,parameters[1],np.double(F0),type(parameters[1])
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
            F3_lst = stdpar[i].split(' ');F4 = [x for x in F3_lst if x is not ''][1]
            parameters[5] = np.longdouble(F4)
    print 'parameters=',parameters,type(parameters[0]),type(parameters[1])
    return parameters

def read_data(datalistname):
    data = np.array([])
    datalist = open(datalistname,'r')
    for file in datalist:
        file = file[0:-1]
        print file
        hdulist = pf.open(file)
        try:
            tdb = hdulist[1].data.field('tdb')
        except:
            print "NO TDB column"
            continue
        data = np.append(data, tdb)
        print len(data)
    data.sort()
    return data

def pfold(time,parfile,duration,f0_flag=True,f1_flag=True,f2_flag=True,f3_flag=False,f4_flag=False,fig_flag=True,bin_profile=1000,threshold=1e6):
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

    # standard profile
    std_pro = np.loadtxt('he_pro_std_1000.dat')
    std_pro = [(x - min(std_pro))/(max(std_pro)-min(std_pro)) for x in std_pro] # Normalization

    # seperate data
    edges = np.arange(min(time),max(time),duration)
    p_x = []
    profile = []
    toa = []
    for i in xrange(len(edges[0:-1])):
        edge0 = edges[i]
        edge1 = edges[i+1]
        data = time[ (time>=edge0) & (time<=edge1) ]
        if len(data)==0:
            print "EMPTY"
            continue
        if len(data)<=threshold:continue
        t0 = min(data)
        T0 = t0/86400 + MJDREFF + MJDREFI
        dt = t0 - pepoch 
        f0 = F0 + F1*dt + (1/2)*F2*(dt**2) + (1/6)*F3*(dt**3) + (1/24)*F4*(dt**4)
        f1 = F1 + F2*dt + (1/2)*F3*(dt**2) + (1/6)*F4*(dt**4)
        f2 = F2 + F3*dt + (1/2)*F4*(dt**2)
        f3 = F3 + F4*dt
        f4 = F4
        print f0,f1,f2

        phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4,1.0)
        p_num = np.histogram(phi,bin_profile)[0]
        p_num = [(x - min(p_num))/(max(p_num)-min(p_num)) for x in p_num] # Normalization
        p_num_x = np.arange(0.,bin_profile,1)/bin_profile
        print len(p_num)

        #toa calculation
        ## ccf shift
        y, delay = ccf(std_pro,p_num)
        y, delay = ccf(p_num,std_pro)
        p_num_std = np.roll(std_pro,delay)
        p_num_std_2 = np.append(p_num_std,p_num_std)
        #

        #p_num_x_2_tmp = p_num_x + 1;p_num_x_2_tmp.tolist();
        #p_num_x_2 = p_num_x.tolist();p_num_2 = p_num
        #p_num_x_2.extend(p_num_x_2_tmp);p_num_2.extend(p_num_2);
        #std_pro_2 = np.append(std_pro,std_pro)
        p_num_x_2_tmp = p_num_x + 1
        p_num_x_2 = np.append(p_num_x,p_num_x_2_tmp)
        p_num_2 = np.append(p_num,p_num);
        std_pro_2 = np.append(std_pro,std_pro)

        if fig_flag:
            plt.figure()
            plt.plot(p_num_x_2,p_num_2,'b')
            plt.plot(p_num_x_2,p_num_std_2,'r')
        
        #phi_peak = p_num_x[p_num.index(max(p_num))]
        phi_peak = p_num_x[np.where(p_num_std==max(p_num_std))][0]
        toa_tmp = T0 + (1/f0) * phi_peak

        p_x.append(p_num_x_2)
        profile.append(p_num_2)
        toa.append(toa_tmp)

    return p_x,profile,toa

def fsearch(time,parfile,duration,fstep,frange,f0_flag=True,f1_flag=True,f2_flag=True,f3_flag=False,f4_flag=False,fig_flag=False,bin_cs=20,bin_profile=1000,threshold=3e6):
    MJDREFF = 0.0007660185
    MJDREFI = 55927
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
        f0 = F0 + F1*dt + (1/2)*F2*(dt**2) + (1/6)*F3*(dt**3) + (1/24)*F4*(dt**4)
        
        f1 = F1 + F2*dt + (1/2)*F3*(dt**2) + (1/6)*F4*(dt**4)
        f2 = F2 + F3*dt + (1/2)*F4*(dt**2)
        f3 = F3 + F4*dt
        f4 = F4
        print 'f0,f1,f2=',
        print '{0:.15f}'.format(f0),
        print '{0:.15f}'.format(f1),f2,type(F0),type(f0)

        f = np.arange(f0-frange,f0+frange,fstep)
        print type(f[0])
        chi_square = np.zeros(len(f))
        N = len(data)
        b = N/bin_cs
        for j in range(0,len(f)):
            #phi_tmp = np.mod((data-t0)*f[j] + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4,1.0)
            phi_tmp = np.mod((data-t0)*f[j] + (1/2)*((data-t0)**2)*f1 + ((data-t0)**3)*f2 ,1.0)
            p_num = np.histogram(phi_tmp,bin_cs)[0]
            #bb = b * np.ones(bin_cs)
            #chi_square[j] = np.sum((p_num-bb)**2)/b
            chi_square[j] = np.std(p_num)**2/np.mean(p_num)
            

            percent = float(j)*100/len(f)
            sys.stdout.write(" fsearch complete: %.2f"%percent);
            sys.stdout.write("%\r");
            sys.stdout.flush()
        print chi_square
        fbest = f[np.where(chi_square==max(chi_square))][0]
        fre[i] = fbest
        tf[i] = T0
        fre_err[i] = fstep/(max(chi_square)/bin_cs)
        if fig_flag:
            plt.figure()
            plt.plot(f,chi_square)
    tf = tf[tf>0]
    fre =fre[fre>0]
    fre_err = fre_err[fre_err>0]

#    return p_num_x_2,p_num_2,f,chi_square
    return tf,fre,fre_err

def toa_write(filename,toa):
    with open(filename,'w')as f:
        f.write('FORMAT 1\n')
        for i in xrange(len(toa)):
            wrt_str = 'HE.toa 3000 %.20f 33 bat\n'%(toa[i])
            f.write(wrt_str)

def fre_write(filename,mjd, fre, fre_err):
    with open(filename,'w')as f:
        for i in xrange(len(fre)):
            wrt_str = '%.16f %.16f %.16f \n'%(mjd[i],fre[i],fre_err[i])
            f.write(wrt_str)

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
    y = [ np.sum(delta_f1 * np.roll(delta_f2,x))/(sigma_f1 * sigma_f2) for x in f2range ]
    delay = np.where(y==max(y))[0] 
    return y,delay

