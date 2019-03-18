#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to ME data with HXMTSoftware(V2) 
#################################

from __future__ import division
import argparse
import os
import glob
import numpy as np

# find and create data dir list
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Example: python me_pipeline.py -i /DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --smfovdet --blinddet --hxbary -r 83.63322083 -d 22.014461')
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-I","--inputlist",help="data archived path in list",type=str)
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("-O","--outputlist",help="products archived path in list",type=str)
#parser.add_argument("--detlist",action="store",help="detector list")
#parser.add_argument("--bigfovdet",action="store_true",help="select big fov detectors of three boxes")
#parser.add_argument("--smfovdet",action="store_true",help="select small fov detectors of three boxes")
#parser.add_argument("--blinddet",action="store_true",help="select all blind detectors")
parser.add_argument("--hxbary",action="store_true",help="carry out Barycentric correction")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()


def main():
    aux_dir = product_dir + "/AUX/" # AUX path
    acs_dir = product_dir + "/ACS/" # ACS path
    me_dir = product_dir + "/ME/"   # ME  path
    clean_dir = product_dir + "/ME/cleaned/"  # ME cleaned data path
    tmp_dir = product_dir + "/ME/tmp/" # ME temporary data
    spec_dir = product_dir +"/ME/spectra/" # spectra results path
    lc_dir = product_dir +"/ME/lightcurve/" #light curve results path


    #make directory for data structure
    if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
    if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
    if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
    if not os.path.isdir(me_dir):os.system('mkdir -p '+me_dir)
    if not os.path.isdir(clean_dir):os.system('mkdir -p '+clean_dir)
    if not os.path.isdir(tmp_dir):os.system('mkdir -p '+tmp_dir)
    if not os.path.isdir(spec_dir):os.system('mkdir -p '+spec_dir)
    if not os.path.isdir(lc_dir):os.system('mkdir -p '+lc_dir)
    
    #read filenames
    print data_dir
    filename = sorted(glob.glob(data_dir + '/ME/*ME-Evt_FFFFFF_V[1-9]*'))[-1]
    orbitname = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
    preciseorbitname = orbitname
    attname = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
    tempfilename = sorted(glob.glob(data_dir + '/ME/HXMT*TH*V[1-9]*'))[-1]
    ehkfilename = sorted(glob.glob(data_dir + '/AUX/*EHK*_V[1-9]*'))[-1]
    
    # select good time intervals utilizing HXMT software
    ## pi calculation
    mepical_text = 'mepical evtfile='+filename+' tempfile='+tempfilename+' outfile='+tmp_dir+'me_pi.fits '+\
            ' clobber=yes history=yes'
    print mepical_text
    os.system(mepical_text)
    
    ## gti selection
    megtigen_text = 'megtigen tempfile='+tempfilename+' ehkfile='+ehkfilename+' outfile='+tmp_dir+'me_gti.fits'+\
            ' defaultexpr=NONE expr="ELV>10&&COR>8&&T_SAA>300&&TN_SAA>300&&ANG_DIST<=0.04"'+\
            ' clobber=yes history=yes'
    print megtigen_text
    os.system(megtigen_text) 
    
    ## ME events reconstruction & deadtime calculation
    megrade_text = 'megrade evtfile='+tmp_dir+'me_pi.fits'+' outfile='+tmp_dir+'me_grade.fits deadfile='+tmp_dir+'me_dtime.fits binsize=1'+\
            ' clobber=yes history=yes'
    print megrade_text
    os.system(megrade_text)

    ## new gti selection
    menewgti_text = 'megti '+tmp_dir+'me_grade.fits '+tmp_dir+'me_gti.fits '+tmp_dir+'me_gti.fits'
    print(menewgti_text)
    try:
        os.system(menewgti_text)
    except:
        print "ERROR: couldn't find megti program"
    
    ## select detectors
    det = '0-5,7,12-23,25,30-41,43,48-53 '
    print 'small FoV detector list:',det
    mescreen_text = 'mescreen evtfile='+tmp_dir+'me_grade.fits gtifile='+tmp_dir+'me_gti.fits outfile='+clean_dir+'me_screen_smfov.fits '+\
    ' userdetid="'+det+'" starttime=0 stoptime=0 minPI=0 maxPI=1024 clobber=yes'
    print mescreen_text
    os.system(mescreen_text)
    # spectra generating
    spec_text = 'mespecgen evtfile="'+clean_dir+'me_screen_smfov.fits" outfile="'+\
            spec_dir+'me_spec_smfov" deadfile="'+tmp_dir+'me_dtime.fits" userdetid="'+\
            det+'" starttime=0 stoptime=0 minPI=0 maxPI=1024' 
    print spec_text
    os.system(spec_text)
    
    det = '10,28,46'
    print 'blind detector list:',det
    mescreen_text = 'mescreen evtfile='+tmp_dir+'me_grade.fits gtifile='+tmp_dir+'me_gti.fits outfile='+clean_dir+'me_screen_blind.fits '+\
    ' baddetfile=$HEADAS/refdata/medetectorstatus.fits'+\
    ' userdetid="'+det+'" starttime=0 stoptime=0 minPI=0 maxPI=1024 clobber=yes'
    print mescreen_text
    os.system(mescreen_text)
    # spectra generating
    spec_text = 'mespecgen evtfile="'+clean_dir+'me_screen_blind.fits" outfile="'+\
            spec_dir+'me_spec_blind" deadfile="'+tmp_dir+'me_dtime.fits" userdetid="'+\
            det+'" starttime=0 stoptime=0 minPI=0 maxPI=1023' 
    print spec_text
    os.system(spec_text)
    
    # carry out barycentering correction
    if args.hxbary:
        # carry out hxbary
        ra = args.ra
        dec = args.dec
        hxbary_text = 'hxbary' + ' ' + clean_dir + 'me_screen_blind.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
        print hxbary_text
        os.system(hxbary_text)


    # carry out barycentering correction
    if args.hxbary:
        # carry out hxbary
        ra = args.ra
        dec = args.dec
        hxbary_text = 'hxbary' + ' ' + clean_dir + 'me_screen_smfov.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
        print hxbary_text
        os.system(hxbary_text)
    mebkgmap(product_dir, clean_dir+'me_screen_blind.fits', ehkfilename, tmp_dir+'me_gti.fits', tmp_dir+'me_dtime.fits',
            tempfilename)
    melcgen(lc_dir, clean_dir+'me_screen_smfov.fits', clean_dir+'me_screen_blind.fits', tmp_dir+'me_dtime.fits', ehkfilename, tmp_dir+'me_gti.fits',
            tempfilename, binsize=1, minPI=100, maxPI=600)

def melcgen(lc_dir, screenfile, blindfile, deadfile, ehkfile, gtifile, tempfile, binsize=1, minPI=100, maxPI=600):
    melc_text = 'melcgen evtfile="'+screenfile+'" outfile="'+lc_dir+'me_lc_smfov'+'" deadfile="'+deadfile+\
            '" userdetid="0-5,7,12-23,25,30-41,43,48-53" binsize='+str(binsize)+' starttime=0 stoptime=0'+\
            ' minPI='+str(minPI)+' maxPI='+str(maxPI)+' deadcorr=no clobber=yes'
    print melc_text
    os.system(melc_text)

    listfile = os.path.join(lc_dir,'me_lc.txt')
    create_listfile_text = "ls %s | sort -V > %s"%(os.path.join(product_dir,'ME','lightcurve','me_lc_smfov_g*'),listfile)
    print create_listfile_text
    os.system(create_listfile_text)
    lcbkgmap_text = 'python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/mebkgmap.py lc '+blindfile +' '+ \
            ehkfile + ' ' + gtifile + ' ' +  deadfile + ' ' + tempfile +' '+\
            listfile +' '+str(minPI)+' '+str(maxPI)+' '+os.path.join(lc_dir,'me_lc_bkg')
    print lcbkgmap_text
    os.system(lcbkgmap_text)

def mebkgmap(product_path, blindfile, ehkfile, gtifile, deadfile, tempfile):
    listfile = os.path.join(product_dir,'ME','spectra','me_spec.txt')
    create_listfile_text = "ls %s | sort -V > %s"%(os.path.join(product_dir,'ME','spectra','me_spec_smfov_g*'),listfile)
    print create_listfile_text
    os.system(create_listfile_text)
    specbkgmap_text = 'python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/mebkgmap.py spec '+blindfile +' '+ \
            ehkfile + ' ' + gtifile + ' ' +  deadfile + ' ' + tempfile +' '+\
            listfile + ' 0 1024 '+ os.path.join(product_path,'ME','spectra','me_bkg_spec')
    print specbkgmap_text
    os.system(specbkgmap_text)

if args.inputlist:
    inputfile = open(args.inputlist)
    outputfile= open(args.outputlist)
    for data_dir,product_dir in zip(inputfile,outputfile):
        data_dir = data_dir[0:-1]
        product_dir = product_dir[0:-1]
        main()
elif args.input == None:
    print 'WARNING: no inputs. "python me_pipeline.py -h" see help'
else:
    data_dir = args.input
    product_dir = args.output
    main()

