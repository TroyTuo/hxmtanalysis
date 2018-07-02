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
parser.add_argument("-f","--file",help="input file name")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("--detlist",action="store",help="detector list")
parser.add_argument("--bigfovdet",action="store_true",help="select big fov detectors of three boxes")
parser.add_argument("--smfovdet",action="store_true",help="select small fov detectors of three boxes")
parser.add_argument("--blinddet",action="store_true",help="select all blind detectors")
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
filename = args.file
print filename
aux_dir = product_dir # AUX path
acs_dir = product_dir # ACS path
me_dir = product_dir  # ME  path

#make direction for data structure
if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(me_dir):os.system('mkdir -p '+me_dir)

#read filenames
filename = filename
print filename 
prefix = filename[0:filename.find(".fits")]
print prefix
orbitname    = sorted(glob.glob(data_dir + '/ID345/orbit*'))[-1]
attname      = sorted(glob.glob(data_dir + '/ID345/ID345*'))[-1]
tempfilename = sorted(glob.glob(data_dir + '/ID321/ID321*'))[-1]
gainfilename = commands.getstatusoutput('ls $CALDB/data/hxmt/me/bcf/*gain*')[1]


# generate ehk file utilizing HXMT software
def ehkgen(infilename,outfile_dir):
    prefix = infilename[0:infilename.find(".fits")]
    orbfile = sorted(glob.glob(data_dir + '/ID345/orbit*'))[-1]
    attfile = sorted(glob.glob(data_dir + '/ID345/ID345*'))[-1]
    outfile = outfile_dir+prefix+"_EHK.fits"
    saafile="/home/hxmt/zhaohaish/soft/hxmtehkgen/SAA/SAA.fits"
    text = "hxmtehkgen orbfile="+orbfile+" attfile="+attfile+" outfile="+outfile+\
            " saafile="+saafile+" step_sec=1 clobber=yes"
    print text
    os.system(text)
    return outfile

ehkfilename = ehkgen(filename,product_dir)

# select good time intervals utilizing HXMT software
## pi calculation
mepical_text = 'mepical evtfile='+data_dir+'/ID315/'+filename+' tempfile='+tempfilename+\
        ' outfile='+me_dir+prefix+'_me_pi.fits gainfile='+gainfilename+\
        ' clobber=yes history=yes'
print mepical_text
#os.system(mepical_text)

## gti selection
megtigen_text = 'megtigen tempfile='+tempfilename+' ehkfile='+ehkfilename+\
        ' outfile='+me_dir+prefix+'_me_gti.fits'+\
        ' defaultexpr=NONE expr="ELV>8&&COR>5&&T_SAA>100&&TN_SAA>100"'+\
        ' clobber=yes history=yes'
print megtigen_text
#os.system(megtigen_text) 

## ME events reconstruction & deadtime calculation
megrade_text = 'megrade evtfile='+me_dir+prefix+'_me_pi.fits'+\
        ' outfile='+me_dir+prefix+'_me_grade.fits deadfile='+me_dir+prefix+'_me_dtime.fits binsize=1'+\
        ' clobber=yes history=yes'
print megrade_text
#os.system(megrade_text)

## select detectors
det = ''
if args.bigfovdet:
    print 'function not available'
if args.smfovdet:
    if det == '':
        det = '0-5,7,12-23,25,30-41,43,48-53'
    else:
        det = det + ',0-5,7,12-23,25,30-41,43,48-53'
if args.blinddet:
    if det == '':
        det = '10,28,46'
    else:
        det = det + ',10,28,46'
if args.detlist:
    det = det + args.detlist;
print 'detector_list:',det

mescreen_text = 'mescreen evtfile='+me_dir+prefix+'_me_grade.fits gtifile='+me_dir+prefix+'_me_gti.fits'+\
        ' outfile='+me_dir+prefix+'_me_screen.fits baddetfile=""'+\
        ' userdetid="'+det+'" starttime=0 stoptime=0 minPI=0 maxPI=1024'
print mescreen_text
#os.system(mescreen_text)

## generate spectrum 
mespecgen_text = 'mespecgen evtfile="'+me_dir+prefix+'_me_screen.fits" deadfile="'+me_dir+prefix+'_me_dtime.fits"'+\
        ' userdetid="'+det+'" eventtype=1 starttime=0 stoptime=0 minPI=0 maxPI=1024 clobber=yes'
print mespecgen_text
#os.system(mespecgen_text)
