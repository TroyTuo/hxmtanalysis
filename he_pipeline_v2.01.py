#!/usr/bin/env python
##################################
#Notice: The Pipeline applies to HE data with HXMTSoftware(V2)
#Usage:
#   Parameters:
#       -i/--input, absolute PATH of archived input data
#       -o/--output, absolute PATH of archived output data
#       --hxbary, barycentering correction flag
#       -r/--ra, right ascension of point source 
#       -d/--dec, declination of point source
#   Example:
#       python he_pipeline.py -i /HXMT_DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --hxbary -r 83.6331 -d 22.1446
#################################


import argparse
import os
import glob

# find and create data dir list

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Example: python he_pipeline.py -i /DATA_PATH/ObID/ -o /OUTPUT_PATH/ObID/ --detlist "0-15,17" --blinddet --nai --hxbary -r 83.63322083 -d 22.014461')
parser.add_argument("-i","--input",help="data archived path")
parser.add_argument("-I","--inputlist",help="data archived path in list",type=str)
parser.add_argument("-o","--output",help="products archived path")
parser.add_argument("-O","--outputlist",help="products archived path in list",type=str)
#parser.add_argument("--nai",action="store_true",help="flag for selecting only NaI Events")
#parser.add_argument("--detlist",action="store",help="detector list")
#parser.add_argument("--blinddet",action="store_true",help="select  blind detectors")
parser.add_argument("--hxbary",action="store_true",help="carry out Barycentric correction")
parser.add_argument("-r","--ra",help="right ascension of barycentering correction",type=float)
parser.add_argument("-d","--dec",help="declination of barycentering correction",type=float)
args = parser.parse_args()


def main():
    aux_dir = product_dir + "/AUX/" # AUX path
    acs_dir = product_dir + "/ACS/" # ACS path
    he_dir = product_dir + "/HE/"   # HE  path
    clean_dir = product_dir + "/HE/cleaned/"  # HE cleaned data path
    tmp_dir = product_dir + "/HE/tmp/" # HE temporary data
    spec_dir = product_dir +"/HE/spectra/" # spectra results path
    lc_dir = product_dir + "/HE/lightcurve/" # light curve results path
    rsp_dir = product_dir + "/HE/rsp/" # RSP results path
    
    #make direction for data structure
    if not os.path.isdir(product_dir):os.system('mkdir -p '+product_dir)
    if not os.path.isdir(aux_dir):os.system('mkdir -p ' +aux_dir)
    if not os.path.isdir(acs_dir):os.system('mkdir -p ' +acs_dir)
    if not os.path.isdir(he_dir):os.system('mkdir -p '+he_dir)
    if not os.path.isdir(clean_dir):os.system('mkdir -p '+clean_dir)
    if not os.path.isdir(tmp_dir):os.system('mkdir -p '+tmp_dir)
    if not os.path.isdir(spec_dir):os.system('mkdir -p '+spec_dir)
    if not os.path.isdir(lc_dir):os.system('mkdir -p '+lc_dir)
    if not os.path.isdir(rsp_dir):os.system('mkdir -p '+rsp_dir)
    
    #read filenames
    filename     = sorted(glob.glob(data_dir + '/HE/*HE-Evt_FFFFFF_V[1-9]*'))[-1]
    print(filename)
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
        print("ERROR: couldn't find prcise orbit file")
    # use EHK file in data file
    try:
        ehkfilename = sorted(glob.glob(data_dir + '/AUX/*EHK*_V[1-9]*'))[-1]
    except:
        print("ERROR: couldn't find EHK file")
    
    # select good time intervals utilizing HXMT software
    ## pi calculation
    hepical_text = 'hepical evtfile='+filename+' outfile='+tmp_dir+'he_pi.fits clobber=yes'
    print(hepical_text)
    os.system(hepical_text)
    ## gti selection
    hegtigen_text = 'hegtigen hvfile='+hvfilename+\
    ' tempfile='+tempfilename+' pmfile='+pmfilename+\
    ' outfile='+tmp_dir+'he_gti.fits ehkfile='+ehkfilename+\
    ' defaultexpr=NONE expr="ELV>10&&COR>8&&TN_SAA>300&&T_SAA>300&&ANG_DIST<=0.04"'+\
    ' pmexpr="" clobber=yes history=yes'
    print(hegtigen_text)
    os.system(hegtigen_text) 
    
    ## select good Events
    det = "0-17"
    hescreen_text = 'hescreen evtfile="'+tmp_dir+'he_pi.fits" gtifile="'+tmp_dir+'he_gti.fits" outfile="'+clean_dir+'he_screen.fits"'+\
    ' userdetid="'+det+'" eventtype=1 anticoincidence=yes'+\
    ' starttime=0 stoptime=0 minPI=0 maxPI=255'+\
    ' clobber=yes history=yes'
    print(hescreen_text)
    os.system(hescreen_text)
    # spectra generating
    spec_text = 'hespecgen evtfile="'+clean_dir+'he_screen.fits" outfile="'+\
            spec_dir+'he_spec" deadfile="'+deadfilename+'" userdetid="'+\
            '0;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17" eventtype=1 starttime=0 '+\
            'stoptime=0 minPI=0 maxPI=255'
    print(spec_text)
    os.system(spec_text)

    ## generate RSP file
    herspgen(product_dir, attname, -1, -91)

    # carry out barycentering correction
    if args.hxbary:
        try:
            os.system('cp ' + preciseorbitname + ' ' + acs_dir)
            # carry out hxbary
            ra = args.ra
            dec = args.dec
            hxbary_text = 'hxbary' + ' ' + clean_dir + 'he_screen.fits' + ' ' + preciseorbitname + ' ' + str(ra) + ' ' + str(dec) + ' ' + '2'
            print(hxbary_text)
            os.system(hxbary_text)
        except:
            print('WARNING: NO Precise Orbit file('+data_dir+')')
    
    hespecgen(product_dir, clean_dir+'he_screen.fits', ehkfilename, tmp_dir+'he_gti.fits', deadfilename)
    helcgen(lc_dir, clean_dir+'he_screen.fits', clean_dir+'he_screen.fits', deadfilename, ehkfilename, tmp_dir+'he_gti.fits',
            binsize=1, minPI=25, maxPI=100)

    return


def helcgen(lc_dir, screenfile, blindfile, deadfile, ehkfile, gtifile, binsize=1, minPI=25, maxPI=100):
    helc_text = 'helcgen evtfile="'+screenfile+'" outfile="'+lc_dir+'he_lc'+'" deadfile="'+deadfile+\
            '" userdetid="0-15,17" binsize='+str(binsize)+' starttime=0 stoptime=0'+\
            ' minPI='+str(minPI)+' maxPI='+str(maxPI)+' deadcorr=no clobber=yes'
    print(helc_text)
    os.system(helc_text)

    listfile = os.path.join(product_dir,'HE','lightcurve','he_lc.txt')
    create_listfile_text = "ls %s | sort -V > %s"%(os.path.join(product_dir,'HE','lightcurve','he_lc_g*'),listfile)
    print(create_listfile_text)
    os.system(create_listfile_text)
    lcbkgmap_text = 'python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hebkgmap.py lc %s %s %s %s %s  0 255 %s'%(blindfile, ehkfile, gtifile,
            deadfile, listfile, lc_dir+'he_lc_bkg')
    print(lcbkgmap_text)
    os.system(lcbkgmap_text)


def hespecgen(product_dir, blindfile, ehkfile, gtifile, deadfile):
    listfile = os.path.join(product_dir,'HE','spectra','he_spec.txt')
    create_listfile_text = "ls %s | sort -V > %s"%(os.path.join(product_dir,'HE','spectra','he_spec_g*'),listfile)
    print(create_listfile_text)
    os.system(create_listfile_text)

    specbkgmap_text = 'python /home/hxmt/hxmtsoft2/soft/hxmtsoft-2.01/hxmt/BKG/BldSpec/hebkgmap.py spec '+blindfile +' '+ \
            ehkfile + ' ' + gtifile + ' ' +  deadfile + ' ' + listfile + ' 0 255 '+\
            os.path.join(product_dir,'HE','spectra','he_spec_bkg')
    print(specbkgmap_text)
    os.system(specbkgmap_text)

def herspgen(product_dir, attfile, ra, dec):
    phafilelist = glob.glob(product_dir+'/HE/spectra/he_spec_g*')
    for phafile in phafilelist:
        index = phafile[phafile.find("g"):-4].split("_")[0]
        outfile = os.path.join(product_dir,'HE','rsp','he_rsp_det'+str(index)+'.fits')
        attfile = attfile
        ra = str(ra)
        dec = str(dec)
        rsp_text = "herspgen %s %s %s %s %s clobber=yes"%(phafile, outfile, attfile, ra, dec)
        print(rsp_text)
        os.system(rsp_text)


if args.inputlist:
    inputfile = open(args.inputlist)
    outputfile= open(args.outputlist)
    for data_dir,product_dir in zip(inputfile,outputfile):
        data_dir = data_dir[0:-1]
        product_dir = product_dir[0:-1]
        main()
elif args.input == None:
    print 'WARNING: no inputs. "python he_pipeline.py -h" see help'
else:
    data_dir = args.input
    product_dir = args.output
    main()
