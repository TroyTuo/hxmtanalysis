#!/usr/bin/env python
from function import fsearch
from function import pf
from function import plt
from function import np
from function import ccf
from argparse import ArgumentParser
'''
a pulsar f search program
Usage: python toacal.py --all --cut
--all: for all Events data; --cut: for TOA calculation of time intervals
'''

tstep = 7200
hdulist = pf.open('1509_screen.fits')
tb = hdulist[1].data
screen_data = tb.field(0)
screen_tdb = tb.field(4)
tmax = max(screen_data)
tmin = min(screen_data)
hdulist.close()

time_range = np.arange(tmin,tmax,tstep)
chi_square_all = [0] * len(time_range)
profile_all = [0] * len(time_range)
time_all = [0] * len(time_range)
fbest_all = [0] * len(time_range)
mjd_all = [0] * len(time_range)
toa_all = [0]* len(time_range)

parser = ArgumentParser()
parser.add_argument("-a","--all",action="store_true",help="fsearch for all events[y/n]")
parser.add_argument("-c","--cut",action="store_true",help="fsearch for cutted events[y/n]")
args = parser.parse_args()

#select time
data_new = np.hstack([screen_data,time_range])
data_new.sort()
if args.all:
    plt.figure()
    print "################## step two ####################"
    f0 = 29.6406
    f0 = 6.582
    f1 = -6.7695e-11
    f2 = 1.9587e-21
    fmax = f0 + 1e-3
    fmin = f0 - 1e-4
    fstep = 0.5e-4
    tdb = screen_tdb
    raw_data = screen_data
    t0 = data_new[1]
    data = tdb - t0

    print fmin,fmax
    plt.subplot(2,1,1)
    p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,0,fstep,errorbar=False,fig=False,pannel=False,bin_cs=20,bin_profile=10)
    fbest = f[chi_square.index(max(chi_square))]
    profile_tot = p
    toa = min(raw_data) + ((1/fbest) * p_x[p.index(max(p))])
    toa = t0 + ((1/fbest) * p_x[p.index(max(p))])
    print 'toa: ',toa
    print t0
    plt.plot(p_x,p)
    min_profile_line = [ p_x[x] for x in xrange(len(p)) if p[x]==min(p) ]
    #for i in min_profile_line:
    #    plt.axvline(x=i-0.025,color='r',linewidth=0.5)
    #    plt.axvline(x=i+0.025,color='r',linewidth=0.5)
    plt.subplot(2,1,2)
    plt.plot(f,chi_square)

if args.cut:
    plt.figure()
    for timei in xrange(len(time_range)):
        print "################## data select ####################"
        tindex0=np.where(data_new == time_range[timei])[0][0]
        if time_range[timei] == time_range[-1]:
            tindex1 = -1
        else:
            tindex1=np.where(data_new == time_range[timei+1])[0][0]
        if tindex1-tindex0 <= 50:
            print "EMPTY"
            continue
        tdb = screen_tdb[tindex0:tindex1]

        raw_data = screen_data[tindex0+1:tindex1-1]

        t0 = data_new[tindex0+1]

        #select baried data
        print "################## step one ####################"
        #data_selected = [tdb[x] for x in index]

        MJD = 57946.00529042
        #data = data_selected - data_selected[0]
        data = tdb - t0

        print "################## step two ####################"
        f0 = 6.582
        f1 = -6.795e-11
        f2 = 1.9587e-21
        fmax = f0 + 1e-3
        fmin = f0 - 1e-3
        fstep = 0.5e-4

        print fmin,fmax
    #    plt.figure()
        plt.subplot(2,1,1)
        p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,0,fstep,errorbar=False,fig=False,pannel=False,bin_cs=20,bin_profile=10)
        fbest = f[chi_square.index(max(chi_square))]
        # CCF correction for toa
        profile_temp = p
        delay = max(ccf(profile_tot,profile_temp))[0]
        profile_temp = np.roll(profile_tot,delay)
        ########################
        toa = min(raw_data) + ((1/fbest) * p_x[p.index(max(p))])
        toa = t0 + ((1/fbest) * p_x[p.index(max(p))])
        toa = t0 + ((1/fbest) * p_x[profile_temp.index(max(profile_temp))])
        print 'toa: ',toa
        print t0
        plt.plot(p_x,p)
        min_profile_line = [ p_x[x] for x in xrange(len(p)) if p[x]==min(p) ]
        #for i in min_profile_line:
        #    plt.axvline(x=i-0.025,color='r',linewidth=0.5)
        #    plt.axvline(x=i+0.025,color='r',linewidth=0.5)
        plt.subplot(2,1,2)
        plt.plot(f,chi_square)

        time_all[timei] = min(raw_data)
        profile_all[timei] = p
        chi_square_all[timei] = chi_square
        fbest_all[timei] = fbest
        mjd_all[timei] = time_range[timei]/86400 + 55927
        toa_all[timei] = toa
    print '###############################################'
    print '###############################################'

    print 'profile'
    print profile_all
    print '\n'
    print 'MJD'
    print mjd_all
    print '\n'
    print 'chi_square'
    print chi_square_all
    print '\n'
    print 'TOA'
    print toa_all
    print '###############################################'
    print '###############################################'
    plt.figure()
    profile_all = [ x for x in profile_all if x is not 0 ]
    plt.imshow(profile_all,cmap='hot',interpolation='nearest')
#plt.show()
