#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to HE data with HXMTSoftware(V2)
#Usage:
#   Parameters:
#       -i/--input, absolute PATH of archived input data
#       -o/--output, absolute PATH of archived output data
#       --nai, NaI Events flag(exclude CsI Events)
#       --detlist, detector list
#       --blinddet, flag that generate screen file for blind detector separately
#       --hxbary, barycentering correction flag
#       -r/--ra, right ascension of point source 
#       -d/--dec, declination of point source
#   Example:
#       python he_pipeline.py -i /HXMT_DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --detlist "0-15,17" --blinddet --nai --hxbary -r 83.6331 -d 22.1446
#################################


import argparse
import os
import glob

# find and create data dir list

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("--nai",action="store_true",help="flag for selecting only NaI Events")
parser.add_argument("--detlist",action="store",help="detector list")
parser.add_argument("--blinddet",action="store_true",help="select all blind detectors")
parser.add_argument("--hxbary",action="store_true",help="carry out hxbary and copy Evt file to local directory")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
aux_dir = product_dir + "/AUX/" # AUX path
acs_dir = product_dir + "/ACS/" # ACS path
he_dir = product_dir + "/HE/"   # HE  path

# generate ehk file utilizing HXMT software
def ehkgen(infile_dir,outfile_dir):
    orbfile = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
    attfile = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
    outfile = outfile_dir+"/AUX/EHK.fits"
    leapfile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/leapsec.fits"
    rigidity="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/rigidity_20060421.fits"
    saafile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/SAA/SAA.fits"
    text = "hxmtehkgen orbfile="+orbfile+" attfile="+attfile+" outfile="+outfile+" leapfile="+leapfile+" rigidity="+rigidity+" saafile="+saafile+" step_sec=1 mean_phi=1 mean_theta=1 mean_psi=1 clobber=yes"
    print text
    os.system(text)

#make direction for data structure
if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(he_dir):os.system('mkdir -p '+he_dir)

#read filenames
filename     = sorted(glob.glob(data_dir + '/HE/*HE-Evt_FFFFFF_V[1-9]*'))[-1]
print filename 
orbitname    = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
attname      = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
hvfilename   = sorted(glob.glob(data_dir + '/HE/HXMT*HV_FFFFFF*V[1-9]*'))[-1]
pmfilename   = sorted(glob.glob(data_dir + '/HE/HXMT*PM_FFFFFF*V[1-9]*'))[-1]
deadfilename = sorted(glob.glob(data_dir + '/HE/HXMT*DTime*V[1-9]*'))[-1]
tempfilename = sorted(glob.glob(data_dir + '/HE/HXMT*TH*V[1-9]*'))[-1]

#copy space craft files
cp_orbit_text = 'cp ' + orbitname + ' ' + acs_dir
os.system(cp_orbit_text)
cp_att_text = 'cp ' + attname + ' ' + acs_dir
os.system(cp_att_text)
try:
    preciseorbitname = glob.glob(data_dir + '/ACS/*Orbit*')[0]
    cp_precise_text = 'cp ' + preciseorbitname + ' ' + acs_dir
    os.system(cp_precise_text)
except:
    print "NO prcise orbit file"
# use EHK file in data file
try:
    ehkfilename = sorted(glob.glob(data_dir + '/AUX/*EHK*_V[1-9]*'))[-1]
    cp_ehk_text = 'cp ' + ehkfilename + ' ' + aux_dir+'/EHK.fits'
    print "###### EHK file is copied from database ######"
    os.system(cp_ehk_text)
except:
    ehkgen(data_dir,product_dir)
    ehkfilename = product_dir + '/ACS/EHK.fits'

# select good time intervals utilizing HXMT software
## pi calculation
if not os.path.isfile(he_dir+'he_pi.fits'):
    hepical_text = 'hepical evtfile='+filename+' outfile='+he_dir+'he_pi.fits clobber=yes'
    os.system(hepical_text)
    print hepical_text
## gti selection
hegtigen_text = 'hegtigen hvfile='+hvfilename+\
' tempfile='+tempfilename+' pmfile='+pmfilename+\
' outfile='+he_dir+'he_gti.fits ehkfile='+ehkfilename+\
' defaultexpr=NONE expr="ELV>6&&COR>8&&TN_SAA>100&&T_SAA>100"'+\
' pmexpr="" clobber=yes history=yes'
os.system(hegtigen_text) 
print hegtigen_text

## select good Events
det = ''
if args.detlist:
    det = det + args.detlist
    # select NaI Events
    if args.nai:
        if not os.path.isfile(he_dir+'he_screen_NaI.fits'):
            hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen_NaI.fits"'+\
            ' baddetfile="" userdetid="'+det+'" eventtype=1 anticoincidence=""'+\
            ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
            ' minpulsewidth=54 maxpulsewidth=70'+\
            ' clobber=yes history=yes'
            os.system(hescreen_text)
            print hescreen_text
            # carry out barycentering correction
            if args.hxbary:
                try:
                    os.system('cp ' + preciseorbitname + ' ' + acs_dir)
                    # carry out hxbary
                    ra = args.ra
                    dec = args.dec
                    hxbary_text = 'hxbary' + ' ' + he_dir + 'he_screen_NaI.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
                    print hxbary_text
                    os.system(hxbary_text)
                except:
                    print 'WARNING: NO Precise Orbit file('+data_dir+')'
    else:
        if not os.path.isfile(he_dir+'he_screen.fits'):
            hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen.fits"'+\
            ' baddetfile="" userdetid="'+det+'" eventtype=1 anticoincidence=""'+\
            ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
            ' minpulsewidth=20 maxpulsewidth=70'+\
            ' clobber=yes history=yes'
            os.system(hescreen_text)
            print hescreen_text
            # carry out barycentering correction
            if args.hxbary:
                try:
                    os.system('cp ' + preciseorbitname + ' ' + acs_dir)
                    # carry out hxbary
                    ra = args.ra
                    dec = args.dec
                    hxbary_text = 'hxbary' + ' ' + he_dir + 'he_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
                    print hxbary_text
                    os.system(hxbary_text)
                except:
                    print 'WARNING: NO Precise Orbit file('+data_dir+')'

if args.blinddet:
    det = '16'
    # select NaI Events
    if args.nai:
        if not os.path.isfile(he_dir+'he_screen_NaI_blind.fits'):
            hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen_NaI_blind.fits"'+\
            ' baddetfile="" userdetid="'+det+'" eventtype=1 anticoincidence=""'+\
            ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
            ' minpulsewidth=54 maxpulsewidth=70'+\
            ' clobber=yes history=yes'
            os.system(hescreen_text)
            print hescreen_text
            # carry out barycentering correction
            if args.hxbary:
                try:
                    os.system('cp ' + preciseorbitname + ' ' + acs_dir)
                    # carry out hxbary
                    ra = args.ra
                    dec = args.dec
                    hxbary_text = 'hxbary' + ' ' + he_dir + 'he_screen_NaI_blind.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
                    print hxbary_text
                    os.system(hxbary_text)
                except:
                    print 'WARNING: NO Precise Orbit file('+data_dir+')'
    else:
        if not os.path.isfile(he_dir+'he_screen_blind.fits'):
            hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen_blind.fits"'+\
            ' baddetfile="" userdetid="'+det+'" eventtype=1 anticoincidence=""'+\
            ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
            ' minpulsewidth=20 maxpulsewidth=70'+\
            ' clobber=yes history=yes'
            os.system(hescreen_text)
            print hescreen_text
            # carry out barycentering correction
            if args.hxbary:
                try:
                    os.system('cp ' + preciseorbitname + ' ' + acs_dir)
                    # carry out hxbary
                    ra = args.ra
                    dec = args.dec
                    hxbary_text = 'hxbary' + ' ' + he_dir + 'he_screen_blind.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
                    print hxbary_text
                    os.system(hxbary_text)
                except:
                    print 'WARNING: NO Precise Orbit file('+data_dir+')'
