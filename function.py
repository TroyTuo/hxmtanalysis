import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt
import scipy.signal as ss
import fileinput

# sigma of sinal
def calsigma(data,bk):
    sigma = (data-bk)/np.sqrt(bk)
    return sigma

# select detector
def sel_det(data,det_id,num):
    det_evt = []
    for i in xrange(len(det_id)):
        if det_id[i] == num:
            det_evt.append(data[i])
    return det_evt

# generate lc
def genlc(data,binsize=1,fig=False,rate=True):
    N = (max(data)-min(data))/binsize
    N = int(N)
    lc = np.histogram(data,N)[0]
    if rate:
        lc = lc/binsize # calculate counts rate instead of counts
    lc_time = np.histogram(data,N)[1][0:-1]
    #null = np.where(lc == 0)
    #lc = np.delete(lc,null)
    #lc_time = np.delete(lc_time,null)
    #print null
    if fig:
        plt.figure('light curve')
        plt.plot(lc_time,lc)
        plt.show()
    return lc_time,lc

# generate spec
def genspec(channel):
    spec = np.histogram(channel,255)[0]
    spec_x = np.linspace(0,255,256)[0:-1]
    return spec_x,spec

# Read data
def readdata(filename='/Users/tuoyouli/Desktop/fermi_toa/data/bary_1deg.fits', sat = 'hxmt',pi_flag = False):
    print 'file is ',filename
    hdulist = pf.open(filename)
    tb = hdulist[1].data

    # data format for different instrument (default is for hxmt 1R level data)
    if sat=='hxmt' or sat=='HXMT':
        '''raw_data,det_id,channel,pulse_width,event_type'''
        raw_data = tb.field(0)
        det_id = tb.field(1)
        channel = tb.field(2)
        pulse_width = tb.field(3)
        event_type = tb.field(5)
        if pi_flag:
            pi = tb.field(7)
            return raw_data,det_id,channel,pulse_width,event_type,pi
        hdulist.close()
        return raw_data,det_id,channel,pulse_width,event_type


    if sat=='fermi' or sat=='FERMI':
        raw_data = tb.field(9)
        hdulist.close()
        return raw_data


def fsearch(data,fmin,fmax,f1,f2,fstep,errorbar=False,fig=False,bin_cs=20,bin_profile=20):    
    
    #data = raw_data - min(raw_data);data.sort();
    raw_data = data
    #data = data - data[0]
    t_0 = raw_data[0]
    N = len(data)
    # bin_cs=20 is DOF for chi_square test
    # bin_profile=20 is profile bin size
    b = N/m
    f = np.arange(fmin,fmax,fstep) 
    #f1 = np.arange(f1min,f1max,f1step)
    chi_square = []
    
    
    for i in range(0,len(f)):
        phi_tmp = np.mod(data*f[i] + (data**2)*f1*0.5 + (data**3)*f2/6,1.0)
        p_num = np.histogram(phi_tmp,bin_cs)[0]
        bb = b * np.ones(bin_cs)
        chi_square.append(np.sum((p_num-bb)**2)/b)
    
    
        fbest = f[chi_square.index(max(chi_square))]
        phi = np.mod(data*fbest + + (data**2)*f1*0.5 + (data**3)*f2/6,1.0)
        p_num = np.histogram(phi,bin_profile)[0]
        p_num_x = np.arange(0.,bin_profile,1)/bin_profile
    
        p_num_x_2_tmp = p_num_x + 1;p_num_x_2_tmp.tolist();
        p_num_x_2 = p_num_x.tolist();p_num_2 = p_num.tolist();
        p_num_x_2.extend(p_num_x_2_tmp);p_num_2.extend(p_num_2);
    
    
    #    with open(filename+'_cp_'+str(cut[j])+'_'+str(cut[j+1]),'wt')as file:
    #        for i in range(0,len(chi_square)):
    #            write_str = '%f %f\n'%(f[i],chi_square[i])
    #            file.write(write_str)
    #
    #    with open(filename+'_profile_'+str(cut[j])+'_'+str(cut[j+1]),'wt')as file:
    #        for i in range(0,len(p_num_x_2)):
    #            write_str = '%f %f\n'%(p_num_x_2[i],p_num_2[i])
    #            file.write(write_str)
    
    print "fbest: ",fbest
    print "T: ",1/fbest
   # toa = t_0 + ((1/fbest) * p_num_x[p_num.index(max(p_num))])
    #toa = t_0 + ((1/fbest) * p_num_x[np.where(p_num == max(p_num))]
    print "t_0(the very first arrived photon): ",t_0
    #print 'TOA: ',toa
    print "done"
    
    if errorbar:
        errorp = [] # statistic error for profile counts
        for i in xrange(len(p_num)):
            errorp.append(p_num[i]**0.5)

    if fig:
        if errorbar:
            plt.figure('results')
            plt.subplot(2,1,1)
            plt.title('profile')
            plt.errorbar(p_num_x,p_num,errorp)
            plt.subplot(2,1,2)
            plt.title('chisquare test')
            plt.plot(f,chi_square)
            plt.show()
        else:
            plt.figure('results')
            plt.subplot(2,1,1)
            plt.plot(p_num_x,p_num)
            plt.subplot(2,1,2)
            plt.plot(f,chi_square)
            plt.show()




    return p_num_x_2,p_num_2,f,chi_square
    
    
def tical(data,fig=False,binsize=0.00001): # calculate distribution of time interval
    ti = []
    data.sort()
    for i in xrange(len(data)-1):
        ti.append(data[i+1]-data[i])
    ti_x,ti_y = genlc(ti,binsize=0.00001,rate=False)
    if fig:
        plt.figure('time interval')
        plt.title('time interval')
        plt.xlabel('time interval($\mu s$)')
        plt.ylabel('counts')
        plt.semilogy(1e6*ti_x,ti_y)
    return ti_x,ti_y









