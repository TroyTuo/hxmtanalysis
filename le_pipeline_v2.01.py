#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to LE data with HXMTSoftware(V2) 
#################################


import argparse
import os
import glob
import numpy as np
import commands

# find and create data dir list
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Example: python le_pipeline.py -i /DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --hxbary -r 83.63322083 -d 22.014461')
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-I","--inputlist",help="data archived path in list",type=str)
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("-O","--outputlist",help="products archived path in list",type=str)
parser.add_argument("--hxbary",action="store_true",help="carry out Barycentric correction")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()

def main():
    aux_dir = product_dir + "/AUX/" # AUX path
    acs_dir = product_dir + "/ACS/" # ACS path
    le_dir = product_dir + "/LE/"   # LE  path
    clean_dir = product_dir + "/LE/cleaned/"  # LE cleaned data path
    tmp_dir = product_dir + "/LE/tmp/" # LE temporary data
    spec_dir = product_dir +"/LE/spectra/" # spectra results path
    lc_dir = product_dir +"/LE/lightcurve/" #light curve results path
    rsp_dir = product_dir + "/LE/rsp/" # RSP results path
    
    
    #make direction for data structure
    if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
    if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
    if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
    if not os.path.isdir(le_dir):os.system('mkdir -p '+le_dir)
    if not os.path.isdir(clean_dir):os.system('mkdir -p '+clean_dir)
    if not os.path.isdir(tmp_dir):os.system('mkdir -p '+tmp_dir)
    if not os.path.isdir(spec_dir):os.system('mkdir -p '+spec_dir)
    if not os.path.isdir(lc_dir):os.system('mkdir -p '+lc_dir)
    if not os.path.isdir(rsp_dir):os.system('mkdir -p '+rsp_dir)
    
    #read filenames
    filename = sorted(glob.glob(data_dir + '/LE/*LE-Evt_FFFFFF_V[1-9]*'))[-1]
    instatusfilename = sorted(glob.glob(data_dir + '/LE/*LE-InsStat_FFFFFF_V[1-9]*'))[-1]
    orbitname = sorted(glob.glob(data_dir + '/ACS/*_Orbit_*V[1-9]*'))[-1]
    preciseorbitname = orbitname
    attname = sorted(glob.glob(data_dir + '/ACS/*_Att_*V[1-9]*'))[-1]
    tempfilename = sorted(glob.glob(data_dir + '/LE/HXMT*TH*V[1-9]*'))[-1]
    gainfilename = commands.getstatusoutput('ls $CALDB/data/hxmt/le/bcf/*gain*')[1]
    ehkfilename = sorted(glob.glob(data_dir + '/AUX/*EHK*'))[-1]
    
    # select good time intervals utilizing HXMT software
    ## pi calculation
    lepical_text = 'lepical evtfile='+filename+' tempfile='+tempfilename+' outfile='+tmp_dir+'le_pi.fits'+\
            ' clobber=yes history=yes'
    print lepical_text
    os.system(lepical_text)
    
    ## le reconstruction
    lerecon_text = 'lerecon evtfile='+tmp_dir+'le_pi.fits outfile='+tmp_dir+'le_recon.fits instatusfile='+instatusfilename+' clobber=yes history=yes'
    print lerecon_text
    os.system(lerecon_text)
    
    ## gti selection
    legtigen_text = 'legtigen evtfile='+filename+' instatusfile='+instatusfilename+' tempfile='+tempfilename+' ehkfile='+ehkfilename+\
            ' outfile='+tmp_dir+'le_gti.fits defaultexpr=NONE rangefile=$HEADAS/refdata/lerangefile.fits'+\
            ' expr="ELV>10&&DYE_ELV>40&&COR>8&&T_SAA>=300&&TN_SAA>=300&&ANG_DIST<=0.04"'+\
            ' clobber=yes history=yes'
    print legtigen_text
    os.system(legtigen_text) 
    ## new git selection
    lenewgti_text = 'python /home/hxmt/hxmtsoft2/hxmtsoftv2.01/install/hxmt/x86_64-pc-linux-gnu-libc2.12/HXMTBKG/soft/legti.py '+tmp_dir+'le_recon.fits '+tmp_dir+'le_gti.fits '+tmp_dir+'le_gti.fits'
    print(lenewgti_text)
    try:
        os.system(lenewgti_text)
        pass
    except:
        print "ERROR: couldn't find legti program"
    
    ## detector selection
    det =   '0,2-4,6-10,12,14,20,22-26,28-30,'+\
            '32,34-36,38-42,44,46,52,54-58,60-62,'+\
            '64,66-68,70-74,76,78,84,86-90,92-94,21,53,85,13,45,77'
    print "small FoV detector list:",det
    ## select good event data
    lescreen_text = 'lescreen evtfile='+tmp_dir+'le_recon.fits gtifile='+tmp_dir+'le_gti.fits outfile='+clean_dir+'le_screen.fits userdetid="'+det+'"'+\
            ' eventtype=1 starttime=0 stoptime=0 minPI=0 maxPI=1535'+\
            ' clobber=yes history=yes'
    print lescreen_text
    os.system(lescreen_text)
    # spectra generating
    spec_text = 'lespecgen evtfile="'+clean_dir+'le_screen.fits" outfile="'+\
            spec_dir+'le_spec" eventtype=1 userdetid="'+\
            '0,2-4,6-10,12,14,20,22-26,28-30,'+\
            '32,34-36,38-42,44,46,52,54-58,60-62,'+\
            '64,66-68,70-74,76,78,84,86-90,92-94'+\
            '" starttime=0 stoptime=0 minPI=0 maxPI=1535' 
    print spec_text
    os.system(spec_text)
    
    ## carry out barycentering correction
    if args.hxbary:
        # carry out hxbary
        ra = args.ra
        dec = args.dec
        hxbary_text = 'hxbary' + ' ' + clean_dir + 'le_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
        print hxbary_text
        os.system(hxbary_text)

    ## generate RSP file
    lerspgen(product_dir, attname, tempfilename, -1, -91)
    

    lebkgmap(product_dir, clean_dir + 'le_screen.fits', tmp_dir + 'le_gti.fits')
    lelcgen(lc_dir, clean_dir+'le_screen.fits', clean_dir+'le_screen.fits', ehkfilename, tmp_dir+'le_gti.fits',
            tempfilename, binsize=1, minPI=150, maxPI=850)

def lelcgen(lc_dir, screenfile, blindfile, ehkfile, gtifile, tempfile, binsize=1, minPI=150, maxPI=850):
    det =   '0,2-4,6-10,12,14,20,22-26,28-30,'+\
            '32,34-36,38-42,44,46,52,54-58,60-62,'+\
            '64,66-68,70-74,76,78,84,86-90,92-94;'
    lelc_text = 'lelcgen evtfile="'+screenfile+'" outfile="'+lc_dir+'le_lc'+'" '\
            ' userdetid="'+det+'" binsize='+str(binsize)+' starttime=0 stoptime=0'+\
            ' minPI='+str(minPI)+' maxPI='+str(maxPI)+' eventtype=1 clobber=yes'
    print lelc_text
    os.system(lelc_text)

    listfile = os.path.join(lc_dir,'le_lc.txt')
    create_listfile_text = "ls %s | sort -V > %s"%(os.path.join(lc_dir,'le_lc_g*'),listfile)
    print create_listfile_text
    os.system(create_listfile_text)
    lcbkgmap_text = 'python /home/hxmt/hxmtsoft2/hxmtsoftv2.01/install/hxmt/x86_64-pc-linux-gnu-libc2.12/HXMTBKG/soft/lebkgmap.py lc '+blindfile +' '+ \
            gtifile + ' ' + listfile +' '+str(minPI)+' '+str(maxPI)+' '+ os.path.join(lc_dir, 'le_lc_bkg')
    print lcbkgmap_text
    os.system(lcbkgmap_text)

def lebkgmap(product_path, blindfile, gtifile):
    listfile = os.path.join(product_dir,'LE','spectra','le_spec.txt')
    create_listfile_text = "ls %s | sort -V > %s"%(os.path.join(product_dir,'LE','spectra','le_spec_g*'),listfile)
    print create_listfile_text
    os.system(create_listfile_text)
    specbkgmap_text = 'python /home/hxmt/hxmtsoft2/hxmtsoftv2.01/install/hxmt/x86_64-pc-linux-gnu-libc2.12/HXMTBKG/soft/lebkgmap.py spec '+blindfile +' '+ \
            gtifile + ' ' + listfile + ' 0 1535 '+ os.path.join(product_path,'LE','spectra','le_bkg_spec')
    print specbkgmap_text
    os.system(specbkgmap_text)

def lerspgen(product_path, attfile, tempfile, ra, dec):
    phafile = glob.glob(product_path+"/LE/spectra/le_spec_g*")[0]
    outfile = os.path.join(product_dir,'LE','rsp','le_rsp.fits')
    attfile = attfile
    tempfile = tempfile
    ra = str(ra)
    dec = str(dec)
    rsp_text = "lerspgen %s %s %s %s %s %s clobber=yes"%(phafile, outfile, attfile, tempfile, ra, dec)
    print(rsp_text)
    os.system(rsp_text)

if args.inputlist:
    inputfile = open(args.inputlist)
    outputfile= open(args.outputlist)
    for data_dir,product_dir in zip(inputfile,outputfile):
        data_dir = data_dir[0:-1]
        product_dir = product_dir[0:-1]
        try:
            main()
        except:
            continue
elif args.input == None:
    print 'WARNING: no inputs. "python le_pipeline.py -h" see help'
else:
    data_dir = args.input
    product_dir = args.output
    main()
