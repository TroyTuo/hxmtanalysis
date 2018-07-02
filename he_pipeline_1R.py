#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to HE data with HXMTSoftware(V2)
#Usage:
#   Parameters:
#       -i/--input, absolute PATH of archived input data
#       -o/--output, absolute PATH of archived output data
#       --nai, NaI Events flag(exclude CsI Events)
#       --detlist, detector list
#   Example:
#       python he_pipeline.py -i /HXMT_DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --detlist "0-15,17" --blinddet --nai --hxbary -r 83.6331 -d 22.1446
#################################


import argparse
import os
import glob

# find and create data dir list

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-f","--file",help="input file name")
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("--nai",action="store_true",help="flag for selecting only NaI Events")
parser.add_argument("--detlist",action="store",help="detector list")
args = parser.parse_args()
data_dir = args.input
product_dir = args.output
filename = args.file
print filename
aux_dir = product_dir # AUX path
acs_dir = product_dir # ACS path
he_dir = product_dir  # HE  path

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

#make direction for data structure
if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
if not os.path.isdir(he_dir):os.system('mkdir -p '+he_dir)

#read filenames
filename     = filename
print filename 
prefix = filename[0:filename.find(".fits")]
print prefix
orbitname    = sorted(glob.glob(data_dir + '/ID345/orbit*'))[-1]
attname      = sorted(glob.glob(data_dir + '/ID345/ID345*'))[-1]
hvfilename   = sorted(glob.glob(data_dir + '/ID310/ID310*'))[-1]
pmfilename   = sorted(glob.glob(data_dir + '/ID311/ID311*'))[-1]
deadfilename = sorted(glob.glob(data_dir + '/ID309/ID309*'))[-1]
tempfilename = sorted(glob.glob(data_dir + '/ID307/ID307*'))[-1]

##copy space craft files
#cp_orbit_text = 'cp ' + orbitname + ' ' + acs_dir
#os.system(cp_orbit_text)
#cp_att_text = 'cp ' + attname + ' ' + acs_dir
#os.system(cp_att_text)
#try:
#    preciseorbitname = glob.glob(data_dir + '/ACS/*Orbit*')[0]
#    cp_precise_text = 'cp ' + preciseorbitname + ' ' + acs_dir
#    os.system(cp_precise_text)
#except:
#    print "NO prcise orbit file"
## use EHK file in data file
#try:
#    ehkfilename = sorted(glob.glob(data_dir + '/AUX/*EHK*_V[1-9]*'))[-1]
#    cp_ehk_text = 'cp ' + ehkfilename + ' ' + aux_dir+'/EHK.fits'
#    print "###### EHK file is copied from database ######"
#    os.system(cp_ehk_text)
#except:

##generate EHK file
ehkfilename = ehkgen(filename,product_dir)

# select good time intervals utilizing HXMT software
## pi calculation
hepical_text = 'hepical evtfile='+data_dir+'/ID301/'+filename+' outfile='+he_dir+prefix+'_he_pi.fits'+\
        ' gainfile=/home/hxmt/hxmtsoft2/CALDB/data/hxmt/he/bcf/hxmt_he_gain_20171030_v1.fits clobber=yes'
print hepical_text
os.system(hepical_text)

## gti selection
hegtigen_text = 'hegtigen hvfile='+hvfilename+\
' tempfile='+tempfilename+' pmfile='+pmfilename+\
' outfile='+he_dir+prefix+'_he_gti.fits ehkfile='+ehkfilename+\
' defaultexpr=NONE expr="ELV>6&&COR>8&&TN_SAA>100&&T_SAA>100"'+\
' pmexpr="" clobber=yes history=yes'
print hegtigen_text
os.system(hegtigen_text) 

## select good Events
det = ''
if args.detlist:
    det = det + args.detlist
    # select NaI Events
    if args.nai:
        hescreen_text = 'hescreen evtfile="'+he_dir+prefix+'_he_pi.fits" gtifile="'+he_dir+prefix+'_he_gti.fits" outfile="'+he_dir+prefix+'_he_screen_NaI.fits"'+\
        ' baddetfile="" userdetid="'+det+'" eventtype=1 anticoincidence=""'+\
        ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
        ' minpulsewidth=54 maxpulsewidth=70'+\
        ' clobber=yes history=yes'
        print hescreen_text
        os.system(hescreen_text)
    else:
        hescreen_text = 'hescreen evtfile="'+he_dir+prefix+'_he_pi.fits" gtifile="'+he_dir+prefix+'_he_gti.fits" outfile="'+he_dir+prefix+'_he_screen.fits"'+\
        ' baddetfile="" userdetid="'+det+'" eventtype=1 anticoincidence=""'+\
        ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
        ' minpulsewidth=20 maxpulsewidth=70'+\
        ' clobber=yes history=yes'
        print hescreen_text
        os.system(hescreen_text)

    ## generate spectrum 
    hespecgen_text = 'hespecgen evtfile="'+he_dir+prefix+'_he_screen.fits" deadfile="'+\
            deadfilename+'" userdetid="'+det+'" eventtype=1 starttime=0 stoptime=0 minPI=0 maxPI=255 clobber=yes'
    print hespecgen_text
    os.system(hespecgen_text)
