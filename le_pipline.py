#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to LE data with HXMTSoftware(V1) 
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
parser.add_argument("--rangefile",action="store_true",help="rangefile path of hegtigen is $hxmtanalysis")
parser.add_argument("--detlist",nargs='+', type=int)
parser.add_argument("--bigfovflag",action="store",help="exclude big fov flag")
parser.add_argument("--midfovflag",action="store",help="exclude mid fov flag")
parser.add_argument("--smfovflag",action="store",help="exclude small fov flag")
parser.add_argument("--blindfovflag",action="store",help="exclude blind fov flag")
parser.add_argument("--gradeflag",action="store",help="upper limit of grade")
parser.add_argument("--hxbary",action="store_true",help="carry out hxbary and copy Evt file to local directory")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
aux_dir = product_dir + "/AUX/" # AUX path
acs_dir = product_dir + "/ACS/" # ACS path
le_dir = product_dir + "/LE/"   # HE  path

#make direction for data structure
if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(le_dir):os.system('mkdir -p '+le_dir)

#read filenames
filename = sorted(glob.glob(data_dir + '/LE/*LE-Evt_FFFFFF_V[1-9]*'))[-1]
orbitname = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
preciseorbitname = sorted(glob.glob(data_dir + '/ACS/*Precise*V[1-9]*'))[-1]
attname = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
tempfilename = sorted(glob.glob(data_dir + '/LE/HXMT*TH*V[1-9]*'))[-1]
print tempfilename 
ehkfilename = aux_dir + "EHK.fits"

#copy space craft files
cp_orbit_text = 'cp ' + orbitname + ' ' + acs_dir
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
#ehkgen(data_dir,product_dir)

# eventype selection
print "select eventype"
hdulist = pf.open(filename)
time = hdulist[1].data.field('time')
et = hdulist[1].data.field('event_type')
hdulist.close()
det = [0,1,2]
if args.detlist:det=args.detlist
detnum = len(det)
etindex = np.equal(et,1)
ettime = time[etindex]
etbin = np.arange(time[0],time[-1],1)
et1rate,et1time = np.histogram(ettime,bins=etbin)
et1time = et1time[:-1]
et1time_tmp = et1time[:-1]
evtindex = np.greater_equal(et1rate,detnum*1000)
edgeindex = evtindex[1:]-evtindex[:-1]
edges = et1time_tmp[edgeindex]
if evtindex[0]:edges=np.concatenate([[time[0]],edges])
if evtindex[-1]:edges=np.concatenate([edges,[time[-1]]])
low_edge_tmp = edges[::2]
up_edge_tmp = edges[1::2]
low_edge = [x for x,y in zip(low_edge_tmp,up_edge_tmp) if y-x >= 2]
up_edge = [y for x,y in zip(low_edge_tmp,up_edge_tmp) if y-x >= 2]
f = open(le_dir+'evttype.tim','w')
for i in xrange(len(low_edge)):
    if i==len(low_edge)-1:
        wrt_str = '(time >= %s && time <= %s) '%(str(low_edge[i]),str(up_edge[i]))
    else:
        wrt_str = '(time >= %s && time <= %s) ||'%(str(low_edge[i]),str(up_edge[i]))
    f.write(wrt_str)
f.close()
# select evttype=1 events
fselect_text = 'fselect '+filename+' '+le_dir+'le_evt.fits @'+le_dir+'evttype.tim clobber=yes'
print fselect_text
#os.system(fselect_text)
filename = le_dir + 'le_evt.fits'

# select good time intervals utilizing HXMT software
## pi calculation
lepical_text = 'lepical evtfile='+filename+' tempfile='+tempfilename+' outfile='+le_dir+'le_pi.fits clobber=yes history=yes'
print lepical_text
#os.system(lepical_text)

## gti selection
legtigen_text = 'legtigen ehkfile='+ehkfilename+' outfile='+le_dir+'le_gti.fits'+\
' defaultexpr=no ELV=15 dye_elv=40 COR=8 sun_angle=10 '+\
'moon_angle=10 SAA=yes T_SAA=100 TN_SAA=100 ang_dist=359'
print legtigen_text
#os.system(legtigen_text) 

## detector selection
bigfov_flag = 'yes'
midfov_flag = 'yes'
smfov_flag = 'no'
blindfov_flag = 'yes'
if args.detlist:det = args.detlist
if args.bigfovflag:bigfov_flag = args.bigfovflag
if args.midfovflag:midfov_flag = args.midfovflag
if args.smfovflag:smfov_flag = args.smfovflag
if args.blindfovflag:blindfov_flag = args.blindfovflag
print bigfov_flag
ledetectorsel_text = 'ledetectorsel infile='+filename+' outfile='+le_dir+'ledetconfig'+' detectorsel="'+str(det)[1:-1]+'" isbigfov='+bigfov_flag+' ismidfov='+midfov_flag+\
' issmlfov='+smfov_flag+' isblindfov='+blindfov_flag
print ledetectorsel_text
#os.system(ledetectorsel_text)

## LE events reconstruction
lerecon_text = 'lerecon infile='+le_dir+'le_pi.fits'+' outfile='+le_dir+'le_grade.fits clobber=yes history=yes'
print lerecon_text
#os.system(lerecon_text)

## select good Events
grade_flag=1
if args.gradeflag:grade_flag=args.gradeflag
lescreen_text = 'lescreen infile='+le_dir+'le_grade.fits'+' userconfigfile='+le_dir+'ledetconfig'+\
' badfile="" gtifile='+le_dir+'le_gti.fits outfile='+le_dir+'le_screen.fits tstart=0 tstop=0 gradeselection="Grade<='+str(grade_flag)+'" otherselection="PHA>0" clobber=yes'
print lescreen_text
#os.system(lescreen_text)

# carry out barycentering correction
if args.hxbary:
    try:
        preciseorbitname = sorted(glob.glob(data_dir + '/ACS/*Precise*V[1-9]*'))[-1]
    except:
        print 'WARNING: NO Precise Orbit file('+data_dir+')'
    # change le_screen keywords and headfile
    cphead_text = 'cphead '+filename+'[Events] '+le_dir+'le_screen.fits[Evt]'
    os.system(cphead_text)
    # carry out hxbary
    ra = args.ra
    dec = args.dec
    hxbary_text = 'hxbary' + ' ' + le_dir + 'le_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
    print hxbary_text
    os.system(hxbary_text)

