#!/usr/bin/env python
from function import *
import argparse

# read profile and parameters
MJD = 55927
t0 = MJD*86400
phase_min = 0.5
phase_max = 0.6
fbest = 29.6397
infile = 'Evt_select.fits'
outfile = 'Evt_phase.fits'
phase_off_file = 'phase_off.fits'
phase_on_file = 'phase_on.fits'

parser = argparse.ArgumentParser()
parser.add_argument("-f0","--fbest",help="best frequency",type=float)
parser.add_argument("-i","--infile",help="name of input file")
parser.add_argument("-pmin","--PhaseMin",help="min edge of phase off bin",type=float)
parser.add_argument("-pmax","--PhaseMax",help="max edge of phase off bin",type=float)
parser.add_argument("-t0","--Epoch",help="best frequency",type=float)
args = parser.parse_args()
t0 = args.Epoch
phase_min = args.PhaseMin
phase_max = args.PhaseMax
fbest = args.fbest

# calculate phase for each events
calc_text = 'ftcalc '+infile+' '+outfile+' Phase "((Time-'+str(t0)+ ')*'+str(fbest)+')%1" clobber=yes'
print calc_text
os.system(calc_text)

# select photons with specific phase
phase_select_text = 'fselect '+outfile+' '+phase_off_file+' "Phase >= '+str(phase_min)+' && Phase <= '+str(phase_max)+'" clobber=yes'
print phase_select_text
os.system(phase_select_text)
phase_select_text = 'fselect '+outfile+' '+phase_on_file+' "Phase <= '+str(phase_min)+' || Phase >= '+str(phase_max)+'" clobber=yes'
print phase_select_text
os.system(phase_select_text)


# generate spectrum
