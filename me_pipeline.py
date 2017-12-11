#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to ME data with HXMTSoftware(V1) 
#################################


import argparse
import os
import glob
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt

# find and create data dir list
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("-c","--config",help="path of ME configure file")
parser.add_argument("--gradeflag",action="store",help="upper limit of grade")
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
preciseorbitname = sorted(glob.glob(data_dir + '/ACS/*Precise*V[1-9]*'))[-1]
attname = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
tempfilename = sorted(glob.glob(data_dir + '/ME/HXMT*TH*V[1-9]*'))[-1]
print tempfilename 
ehkfilename = aux_dir + "EHK.fits"

#copy space craft files cp_orbit_text = 'cp ' + orbitname + ' ' + acs_dir
if not glob.glob(acs_dir+'*_Orbit_*V*'):os.system(cp_orbit_text)
cp_att_text = 'cp ' + attname + ' ' + acs_dir
if not glob.glob(acs_dir+'*_Att_*V*'):os.system(cp_att_text)
cp_precise_text = 'cp ' + preciseorbitname + ' ' + acs_dir
if not glob.glob(acs_dir+'*Precise*'):os.system(cp_precise_text)

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
mepical_text = 'mepical evtfile='+filename+' tempfile='+tempfilename+' outfile='+me_dir+'me_pi.fits clobber=yes history=yes'
print mepical_text
os.system(mepical_text)

## gti selection
megtigen_text = 'megtigen ehkfile='+ehkfilename+' outfile='+me_dir+'me_gti.fits'+\
' defaultexpr=no ELV=15 COR=8 sun_angle=10 '+\
'moon_angle=10 SAA=yes T_SAA=10 TN_SAA=10 ang_dist=359'
print megtigen_text
os.system(megtigen_text) 

## ME events reconstruction
megrade_text = 'megrade infile='+me_dir+'me_pi.fits'+' outfile='+me_dir+'me_grade.fits clobber=yes history=yes'
print megrade_text
os.system(megrade_text)


## select good Events
grade_flag=1
configname = args.config
if args.gradeflag:grade_flag=args.gradeflag
mescreen_text = 'mescreen infile='+me_dir+'me_grade.fits gtifile='+me_dir+'me_gti.fits outfile='+me_dir+'me_screen.fits outdetfile='+\
me_dir+'detconfig confile='+configname+' grade="grade<='+str(grade_flag)+'" baddetfile=""'
print mescreen_text
os.system(mescreen_text)

# carry out barycentering correction
if args.hxbary:
    try:
        preciseorbitname = sorted(glob.glob(data_dir + '/ACS/*Precise*V[1-9]*'))[-1]
    except:
        print 'WARNING: NO Precise Orbit file('+data_dir+')'
    # change le_screen keywords and headfile
#    cphead_text = 'cphead '+filename+'[Events] '+me_dir+'me_screen.fits[Evt]'
#    os.system(cphead_text)
    # carry out hxbary
    ra = args.ra
    dec = args.dec
    hxbary_text = 'hxbary' + ' ' + me_dir + 'me_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
    print hxbary_text
    os.system(hxbary_text)

