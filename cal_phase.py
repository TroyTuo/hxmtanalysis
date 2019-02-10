#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
from astropy.io import fits
from fitsio import FITS
import numpy as np 
import matplotlib.pyplot as plt
import sys
import argparse


#read data file
def read_data(filename, colname='TDB'):
    hdulist = fits.open(filename)
    tdb = hdulist[1].data.field(colname)
    hdulist.close()
    return tdb

#read par file
def read_par(parname):
    pardata = open(parname,'r')
    stdpar = []
    parameters = np.zeros(13,dtype='longdouble')
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

def check_subset(time, tstart, tfinish):
    if min(time) >= tstart and tfinish:
        pass
        #print("time is a subset of [START, FINISH]")
    else:
        print("WARNING: time set is not in parameter time range [START, FINISH]")

def phi_cal(time, parfile, instrument='hxmt'):
    if instrument == 'hxmt' or instrument == 'HXMT':
        MJDREFF = 0.0007660185
        MJDREFI = 55927
    elif instrument == 'fermi' or instrument == 'FERMI':
        MJDREFF = 0.00074287037037037
        MJDREFI = 51910

    #read parfile and parameters
    parameters = read_par(parfile)
    check_subset(time, parameters[11], parameters[12])
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
    dt = t0 - pepoch 
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

    phi = np.mod((data-t0)*f0 + (1/2)*((data-t0)**2)*f1 + (1/6)*((data-t0)**3)*f2 + (1/24)*((data-t0)**4)*f3 + (1/120)*((data-t0)**5)*f4 +
            (1/np.math.factorial(6))*((data-t0)**6)*f5 + (1/np.math.factorial(7))*((data-t0)**7)*f6 + (1/np.math.factorial(8))*((data-t0)**8)*f7 + 
            (1/np.math.factorial(9))*((data-t0)**9)*f8 + (1/np.math.factorial(10))*((data-t0)**10)*f9 ,1.0)
    #print phi
    print "...processing..."
    return phi

def write_file(datafile, phi):
    try:
        hdulist = FITS(datafile,'rw')
        hdulist[1].insert_column(name='Phase', data=phi)
        hdulist.close()
    except:
        hdulist = fits.open(datafile)
        table = hdulist[1].data
        table['Phase'][:] = phi
        hdulist.writeto(datafile, overwrite=True)
    print "...adding a column to event file..."
    print "Success"
    return

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
            description='Example: python cal_phase.py -f eventfile.FITS -p ephemeris.par')
    parser.add_argument("-f","--evtfile",help="The Event file containing the column of the Barycenter corrected time")
    parser.add_argument("-p","--parfile",help="The ephemeris file")
    parser.add_argument("--colname", help='The column name of Barycenter corrected time(the default value is "TDB"')
    parser.add_argument("--instrument", help='The name of Instrument(HXMT/FERMI) (the default value is "HXMT"')
    args = parser.parse_args()

    evtfile = args.evtfile
    parfile = args.parfile
    if args.colname:
        colname = args.colname
    else:
        colname = 'TDB'
    if args.instrument:
        instrument = args.instrument
    else:
        instrument = 'HXMT'

    time = read_data(evtfile,colname=colname)
    phi = phi_cal(time, parfile, instrument=instrument)
    write_file(evtfile, phi)
    return


if __name__ == "__main__":
    if len(sys.argv) <=1:
        print "RUN 'python cal_phase.py --help' for help"
    else:
        main()
