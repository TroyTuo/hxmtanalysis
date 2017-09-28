#!/usr/bin/env python

import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt
import scipy.signal as ss
import fileinput
import commands
import os
import sys
from time import sleep


# sigma of sinal
def calsigma(data,bk):
    sigma = (data-bk)/np.sqrt(bk)
    return sigma

# select detector
def sel_det(data,det_id,num):
    det_evt = [det_id[x] for x in xrange(len(det_id)) if det_id[x]==num]


   # det_evt = []
   # for i in xrange(len(det_id)):
   #     if det_id[i] == num:
   #         det_evt.append(data[i])
    return det_evt

# generate lc
def genlc(data,binsize=1,fig=False,rate=True,pannel=True):
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
        plt.plot(lc_time,lc)
    if pannel:
        # annotation text
        text = 'Entries = '+str(len(data))+'\n binsize = '+str(binsize)+'\n duration = '+str(max(lc_time)-min(lc_time))
        plt.annotate(text, xy=(1, 1), xytext=(-15, -15), fontsize=10,
        xycoords='axes fraction', textcoords='offset points',
        bbox=dict(facecolor='white', alpha=0.8),
        horizontalalignment='right', verticalalignment='top')

    return lc_time,lc

# generate spec
def genspec(channel):
    spec = np.histogram(channel,255)[0]
    spec_x = np.linspace(0,255,256)[0:-1]
    return spec_x,spec

# Read data
def readdata(filename='/Users/tuoyouli/Desktop/fermi_toa/data/bary_1deg.fits', sat = 'hxmt',bary_flag=False,pi_flag = False,quiet = True):
    """return raw_data,det_id,channel,pulse_width,acd,event_type
    for 1R level data;
    return raw_data, for FERMI data;
    return raw_data,det_id,channel,pulse_width,acd,event_type,pi
    for 1R level data after pi calculation
    """
    if quiet != True:print 'file is ',filename
    hdulist = pf.open(filename)
    tb = hdulist[1].data

    # data format for different instrument (default is for hxmt 1R level data)
    if sat=='hxmt' or sat=='HXMT':
        '''raw_data,det_id,channel,pulse_width,event_type'''
        raw_data = tb.field(0)
        det_id = tb.field(1)
        channel = tb.field(2)
        pulse_width = tb.field(3)
        acd = tb.field(4)
        event_type = tb.field(5)
        if bary_flag:
            raw_data = tb.field(0)
            tdb = tb.field(1)
            det_id = tb.field(2)
            channel = tb.field(3)
            pulse_width = tb.field(4)
            acd = tb.field(5)
            event_type = tb.field(6)
            hdulist.close()
            return raw_data,tdb,det_id,channel,pulse_width,acd,event_type
        if pi_flag:
            pi = tb.field(7)
            hdulist.close()
            return raw_data,det_id,channel,pulse_width,acd,event_type,pi
        hdulist.close()
        return raw_data,det_id,channel,pulse_width,acd,event_type


    if sat=='fermi' or sat=='FERMI':
        raw_data = tb.field(9)
        hdulist.close()
        return raw_data

def readehk(ehkfilename):
    """time,sat_alt,sat_lon,sat_lat,elv,dye_elv,ang_dist,cor,saa_flag,sun_ang,moon_ang"""
    hdulist = pf.open(ehkfilename)
    tb = hdulist[1].data
    time = tb.field(0)
    sat_alt = tb.field(7)
    sat_lon = tb.field(8)
    sat_lat = tb.field(9)
    elv = tb.field(10)
    dye_elv = tb.field(11)
    ang_dist = tb.field(19)
    cor = tb.field(20)
    saa_flag = tb.field(22)
    sun_ang = tb.field(25)
    moon_ang = tb.field(26)
    hdulist.close()
    return time,sat_alt,sat_lon,sat_lat,elv,dye_elv,ang_dist,cor,saa_flag,sun_ang,moon_ang

def readpm(pmfilename,num=0):
    """time,pm_counts"""
    hdulist = pf.open(pmfilename)
    tb = hdulist[1].data
    time = tb.field(0)
    pm_cnt = tb.field(num+1)
    hdulist.close()
    return time,pm_cnt

def readhv(hvfilename,phonum=0):
    """time,phoDet0-17"""
    hdulist = pf.open(hvfilename)
    tb = hdulist[1].data
    time = tb.field(0)
    pho_cnt = tb.field(1)
    hdulist.close()
    return time,pho_cnt



def fsearch(data,fmin,fmax,f1,f2,fstep,errorbar=False,fig=False,pannel=True,bin_cs=20,bin_profile=20):

    #data = raw_data - min(raw_data);data.sort();
    raw_data = data
    #data = data - data[0]
    t_0 = raw_data[0]
    N = len(data)
    # bin_cs=20 is DOF for chi_square test
    # bin_profile=20 is profile bin size
    b = N/bin_cs
    f = np.arange(fmin,fmax,fstep)
    #f1 = np.arange(f1min,f1max,f1step)
    chi_square = [0] * len(f)


    for i in range(0,len(f)):
        phi_tmp = np.mod(data*f[i] + (data**2)*f1*0.5 + (data**3)*f2/6,1.0)
        p_num = np.histogram(phi_tmp,bin_cs)[0]
        bb = b * np.ones(bin_cs)
        chi_square[i] = (np.sum((p_num-bb)**2)/b)

        percent = float(i)*100/len(f)
        sys.stdout.write(" fsearch complete: %.2f"%percent);
        sys.stdout.write("%\r");
        sys.stdout.flush()
    print '\n'

    fbest = f[chi_square.index(max(chi_square))]
    phi = np.mod(data*fbest + (data**2)*f1*0.5 + (data**3)*f2/6,1.0)
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
    print "t_0(the very first arrived photon): ",t_0
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
    
    if pannel:
        # annotation text
        text = 'duration = '+str(max(raw_data)-min(raw_data))+'\n bin = '+str(bin_profile)+'\n fbest = '+str(fbest)
        plt.annotate(text, xy=(1, 1), xytext=(-15, -15), fontsize=10,
        xycoords='axes fraction', textcoords='offset points',
        bbox=dict(facecolor='white', alpha=0.8),
        horizontalalignment='right', verticalalignment='top')

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

def acddel(data,acd):
    #we exclude events once there is acd triggered
    acdsum = [ x.sum() for x in acd]
    index = [x for x in xrange(len(acdsum)) if acdsum[x]>0]
    selected_data = [data[x] for x in index]
    return selected_data

def ehkgen(infile_dir,outfile_dir):
    v2_flag = commands.getoutput('if [ -f '+infile_dir+'/ACS/*_Orbit_*V2* ];then echo 1;else echo 0;fi')                                                      
    if v2_flag == '1':
        orbfile = commands.getoutput('ls ' + infile_dir + '/ACS/*_Orbit_*V2*')
    else:
        v1_flag = commands.getoutput('if [ -f '+infile_dir+'/ACS/*_Orbit_*V1* ];then echo 1;else echo 0;fi')                                                      
        if v1_flag == '1':
            orbfile = commands.getoutput('ls ' + infile_dir + '/ACS/*_Orbit_*V1*')
        else:
            orbfile = commands.getoutput('ls '+infile_dir+'/ACS/*_Orbit_*')

    attfile= commands.getoutput('ls '+infile_dir+'/ACS/H*Att*')
    outfile = outfile_dir+"/AUX/EHK.fits"
    leapfile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/leapsec.fits"
    rigidity="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/rigidity_20060421.fits"
    saafile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/SAA/SAA.fits"
    text = "hxmtehkgen orbfile="+orbfile+" attfile="+attfile+" outfile="+outfile+" leapfile="+leapfile+" rigidity="+rigidity+" saafile="+saafile+" step_sec=1 mean_phi=0.1 mean_theta=0.1 mean_psi=0.1"
    print text
    os.system(text)

def readscreen(filename):
    hdulist = pf.open(filename)
    tb = hdulist[1].data
    time = tb.field(0)
    pi = tb.field(1)
    det_id = tb.field(2)
    hdulist.close()
    return time,pi,det_id

def readlc(filename):
    hdulist = pf.open(filename)
    tb = hdulist[1].data
    time = tb.field(0)
    rate = tb.field(1)
    error = tb.field(2)
    hdulist.close()
    return time,rate,error

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

