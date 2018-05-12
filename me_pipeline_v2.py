#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to ME data with HXMTSoftware(V2) 
#################################

from __future__ import division
import argparse
import os
import glob
import commands
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt

# find and create data dir list
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("--detlist",action="store",help="detector list")
parser.add_argument("--bigfovdet",action="store_true",help="select big fov detectors of three boxes")
parser.add_argument("--smfovdet",action="store_true",help="select small fov detectors of three boxes")
parser.add_argument("--blinddet",action="store_true",help="select all blind detectors")
parser.add_argument("--hxbary",action="store_true",help="carry out hxbary and copy Evt file to local directory")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
aux_dir = product_dir + "/AUX/" # AUX path
acs_dir = product_dir + "/ACS/" # ACS path
me_dir = product_dir + "/ME/"   # HE  path

#make direction for data structure
if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(me_dir):os.system('mkdir -p '+me_dir)

#read filenames
filename = sorted(glob.glob(data_dir + '/ME/*ME-Evt_FFFFFF_V[1-9]*'))[-1]
orbitname = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
try:
    preciseorbitname = sorted(glob.glob(data_dir + '/ACS/*Orbit*V[1-9]*'))[-1]
except:
    print "NO Precise Orbit file"+data_dir
attname = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
tempfilename = sorted(glob.glob(data_dir + '/ME/HXMT*TH*V[1-9]*'))[-1]
ehkfilename = aux_dir + "EHK.fits"
gainfilename = commands.getstatusoutput('ls $CALDB/data/hxmt/me/bcf/*gain*')[1]

#copy space craft files 
cp_orbit_text = 'cp ' + orbitname + ' ' + acs_dir
if not glob.glob(acs_dir+'*_Orbit_*V*'):
    os.system(cp_orbit_text)
    cp_att_text = 'cp ' + attname + ' ' + acs_dir
if not glob.glob(acs_dir+'*_Att_*V*'):
    os.system(cp_att_text)
try:
    cp_precise_text = 'cp ' + preciseorbitname + ' ' + acs_dir
    if not glob.glob(acs_dir+'*Precise*'):
        os.system(cp_precise_text)
except:
    print "NO precise orbit file"+data_dir

# generate ehk file utilizing HXMT software
def ehkgen(infile_dir,outfile_dir):
    orbfile = sorted(glob.glob(infile_dir+'/ACS/*_Orbit_*V[1-9]*'))[-1]
    attfile = sorted(glob.glob(infile_dir+'/ACS/H*Att*V[1-9]*'))[-1]
    outfile = outfile_dir+"/AUX/EHK.fits"
    leapfile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/leapsec.fits"
    rigidity="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/rigidity_20060421.fits"
    saafile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/SAA/SAA.fits"
    text = "hxmtehkgen orbfile="+orbfile+" attfile="+attfile+" outfile="+outfile+" leapfile="+leapfile+" rigidity="+rigidity+" saafile="+saafile+" step_sec=1 mean_phi=0.1 mean_theta=0.1 mean_psi=0.1 clobber=yes"
    print text
    os.system(text)
if not glob.glob(aux_dir+'/EHK.fits'):
    ehkgen(data_dir,product_dir)

# select good time intervals utilizing HXMT software
## pi calculation
mepical_text = 'mepical evtfile='+filename+' tempfile='+tempfilename+' outfile='+me_dir+'me_pi.fits gainfile='+gainfilename+\
        ' clobber=yes history=yes'
print mepical_text
os.system(mepical_text)

## gti selection
megtigen_text = 'megtigen tempfile='+tempfilename+' ehkfile='+ehkfilename+' outfile='+me_dir+'me_gti.fits'+\
        ' defaultexpr=NONE expr="ELV>8&&COR>5&&T_SAA>100&&TN_SAA>100"'+\
        ' clobber=yes history=yes'
print megtigen_text
os.system(megtigen_text) 

## ME events reconstruction & deadtime calculation
megrade_text = 'megrade evtfile='+me_dir+'me_pi.fits'+' outfile='+me_dir+'me_grade.fits deadfile='+me_dir+'me_dtime.fits binsize=1'+\
        ' clobber=yes history=yes'
print megrade_text
os.system(megrade_text)

## select detectors
det = ''
if args.bigfovdet:
    print 'function not available'
if args.smfovdet:
    det = det + '0-5,7,12-23,25,30-41,43,48-53'
if args.blinddet:
    det = det + '10,28,46'
if args.detlist:
    det = det + args.detlist;
print 'detector_list:',det

mescreen_text = 'mescreen evtfile='+me_dir+'me_grade.fits gtifile='+me_dir+'me_gti.fits outfile='+me_dir+'me_screen.fits baddetfile=""'+\
' userdetid="'+det+'" starttime=0 stoptime=0 minPI=0 maxPI=1024'
print mescreen_text
os.system(mescreen_text)

# carry out barycentering correction
if args.hxbary:
    try:
        preciseorbitname = sorted(glob.glob(data_dir + '/ACS/*Orbit*V[1-9]*'))[-1]
    except:
        print 'WARNING: NO Precise Orbit file('+data_dir+')'
    # carry out hxbary
    ra = args.ra
    dec = args.dec
    hxbary_text = 'hxbary' + ' ' + me_dir + 'me_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
    print hxbary_text
    os.system(hxbary_text)

