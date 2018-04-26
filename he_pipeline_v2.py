#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to HE data with HXMTSoftware(V1) 
#Usage:
#   Parameters:
#       -i/--input, absolute PATH of archived input data
#       -o/--output, absolute PATH of archived output data
#       --nai, NaI Events flag(exclude CsI Events)
#       --hxbary, barycentering correction flag
#       -r/--ra, right ascension of point source 
#       -d/--dec, declination of point source
#   Example:
#       python he_pipeline.py -i /HXMT_DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --nai --hxbary -r 83.6331 -d 22.1446
#################################


import argparse
import os
import glob

# find and create data dir list

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("--rangefile",action="store_true",help="rangefile path of hegtigen is $hxmtanalysis")
parser.add_argument("--nai",action="store_true",help="flag for selecting only NaI Events")
parser.add_argument("--hxbary",action="store_true",help="carry out hxbary and copy Evt file to local directory")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
aux_dir = product_dir + "/AUX/" # AUX path
acs_dir = product_dir + "/ACS/" # ACS path
he_dir = product_dir + "/HE/"   # HE  path

#make direction for data structure
if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(he_dir):os.system('mkdir -p '+he_dir)

#read filenames
filename = glob.glob(data_dir + '/HE/*HE-Evt_FFFFFF_V*')[0]
print filename 
print data_dir
print glob.glob(data_dir + '/ACS/*_Orbit_*V*')
orbitname = glob.glob(data_dir + '/ACS/*_Orbit_*V*')[0]
print data_dir
attname = glob.glob(data_dir + '/ACS/*_Att_*V*')[0]
hvfilename = glob.glob(data_dir + '/HE/HXMT*HV_FFFFFF*')[0]
pmfilename = glob.glob(data_dir + '/HE/HXMT*PM_FFFFFF*')[0]
deadfilename = glob.glob(data_dir + '/HE/HXMT*DTime*')[0]
tempfilename = glob.glob(data_dir + '/HE/HXMT*TH*')[0]
ehkfilename = aux_dir + "EHK.fits"

#copy space craft files
cp_orbit_text = 'cp ' + orbitname + ' ' + acs_dir
os.system(cp_orbit_text)
cp_att_text = 'cp ' + attname + ' ' + acs_dir
os.system(cp_att_text)
try:
    preciseorbitname = glob.glob(data_dir + '/ACS/*Precise*')[0]
    cp_precise_text = 'cp ' + preciseorbitname + ' ' + acs_dir
    os.system(cp_precise_text)
except:
    print "NO prcise orbit file"

# generate ehk file utilizing HXMT software
def ehkgen(infile_dir,outfile_dir):
    orbfile = glob.glob(infile_dir+'/ACS/*_Orbit_*')[0]
    attfile = glob.glob(infile_dir+'/ACS/H*Att*')[0]
    outfile = outfile_dir+"/AUX/EHK.fits"
    leapfile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/leapsec.fits"
    rigidity="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/refdata/rigidity_20060421.fits"
    saafile="/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/SAA/SAA.fits"
    text = "hxmtehkgen orbfile="+orbfile+" attfile="+attfile+" outfile="+outfile+" leapfile="+leapfile+" rigidity="+rigidity+" saafile="+saafile+" step_sec=1 mean_phi=1 mean_theta=1 mean_psi=1 clobber=yes"
    print text
    os.system(text)
ehkgen(data_dir,product_dir)

# select good time intervals utilizing HXMT software
## pi calculation
hepical_text = 'hepical evtfile='+filename+' outfile='+he_dir+'he_pi.fits clobber=yes'
os.system(hepical_text)
print hepical_text
## gti selection
hegtigen_text = 'hegtigen hvfile='+hvfilename+\
' tempfile='+tempfilename+' outfile='+he_dir+'he_gti.fits ehkfile='+ehkfilename+\
' defaultexpr=no ELV=6 COR=10 SAA=yes T_SAA=100 TN_SAA=100 SUN_ANGLE=10'\
' MOON_ANGLE=10 ANG_DIST=100 clobber=yes'
if args.rangefile:
    hegtigen_text = 'hegtigen hvfile='+hvfilename+\
    ' tempfile='+tempfilename+' outfile='+he_dir+'he_gti.fits ehkfile='+ehkfilename+\
    ' defaultexpr=yes rangefile=/home/hxmt/tuoyl/HXMT_process/hxmtanalysis/HEParamConfig.fits'
os.system(hegtigen_text) 
print hegtigen_text
## select good Events
hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen.fits"'\
' baddetfile="" detid="det_id>=0" eventtype=1 anticoincidence=""'
os.system(hescreen_text)
print hescreen_text

# select NaI Events
if args.nai:
    pwselect_text = 'fselect ' + he_dir + 'he_pi.fits ' + he_dir + 'he_pi_NaI.fits "Pulse_Width >= 54 && Pulse_Width <= 70" clobber=yes'
    os.system(pwselect_text)
    
    hescreen_text = 'hescreen evtfile="'+he_dir+'he_pi_NaI.fits" gtifile="'+he_dir+'he_gti.fits" outfile="'+he_dir+'he_screen_NaI.fits"'\
    ' baddetfile="" detid="det_id>=0" eventtype=1 anticoincidence=""'
    os.system(hescreen_text)
    print hescreen_text
    # carry out barycentering correction
    if args.hxbary:
        try:
            preciseorbitname = glob.glob(data_dir + '/ACS/*Precise*')[0]
            os.system('cp ' + preciseorbitname + ' ' + acs_dir)
            # change he_screen keywords and headfile
            fparkey_text = 'fparkey Events '+he_dir+'he_screen_NaI.fits EXTNAME'
            os.system(fparkey_text)
            cphead_text = 'cphead '+filename+'[Events] '+he_dir+'he_screen_NaI.fits[Events]'
            os.system(cphead_text)
            # carry out hxbary
            ra = args.ra
            dec = args.dec
            hxbary_text = 'hxbary' + ' ' + he_dir + 'he_screen_NaI.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
            print hxbary_text
            os.system(hxbary_text)
        except:
            print 'WARNING: NO Precise Orbit file('+data_dir+')'


# carry out barycentering correction
if args.hxbary:
    try:
        preciseorbitname = glob.glob(data_dir + '/ACS/*Precise*')[0]
        os.system('cp ' + preciseorbitname + ' ' + acs_dir)
        # change he_screen keywords and headfile
        fparkey_text = 'fparkey Events '+he_dir+'he_screen.fits EXTNAME'
        os.system(fparkey_text)
        cphead_text = 'cphead '+filename+'[Events] '+he_dir+'he_screen.fits[Events]'
        os.system(cphead_text)
        # carry out hxbary
        ra = args.ra
        dec = args.dec
        hxbary_text = 'hxbary' + ' ' + he_dir + 'he_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
        print hxbary_text
        os.system(hxbary_text)
    except:
        print 'WARNING: NO Precise Orbit file('+data_dir+')'

