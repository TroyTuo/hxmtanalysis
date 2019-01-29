#!/usr/bin/env python

import os
import argparse
import glob

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description='')
parser.add_argument("-r","--rawdata",help="data archived path")
parser.add_argument("-p","--product",help="data archived path")
args = parser.parse_args()
rawdata_path = args.rawdata
product_path = args.product

def readfile(instru='HE'):
    blindfile = glob.glob(os.path.join(product_path,instru,'cleaned','*blind*'))
    ehkfile = glob.glob(os.path.join(rawdata_path,'AUX','*EHK*'))
    gtifile = glob.glob(os.path.join(product_path,instru,'tmp','*gti*'))
    tempfile = glob.glob(os.path.join(rawdata_path,instru,'*TH*'))
    if instru == 'HE':
        deadfile = glob.glob(os.path.join(rawdata_path,instru,'*DTime*'))
        screenfile = glob.glob(os.path.join(product_path,instru,'cleaned','he_screen_NaI.fits'))
    elif instru == 'ME':
        deadfile = glob.glob(os.path.join(product_path,instru,'tmp','*dtime*'))
        screenfile = glob.glob(os.path.join(product_path,instru,'cleaned','me_screen_smfov.fits'))
    else:
        deadfile = []
        screenfile = glob.glob(os.path.join(product_path,instru,'cleaned','le_screen_smfov.fits'))

    return screenfile, blindfile, ehkfile, gtifile, tempfile, deadfile

def helcgen(evtfile, outpath, deadfile, binsize=1, minPI=20, maxPI=200):
    lc_text = 'helcgen evtfile=%s outfile=%s deadfile=%s userdetid=%s eventtype=1 binsize=%s minPI=%s maxPI=%s clobber=yes'\
            %(evtfile[0], outpath+'he_lc', deadfile[0], '"0-15,17"',\
            str(binsize),str(minPI), str(maxPI))
    print(lc_text)
    os.system(lc_text)
    return

def melcgen(evtfile, outpath, deadfile, binsize=1, minPI=0, maxPI=1024):
    lc_text = 'melcgen evtfile=%s outfile=%s deadfile=%s userdetid=%s binsize=%s minPI=%s maxPI=%s clobber=yes'\
            %(evtfile[0], outpath+'me_lc_smfov', deadfile[0], '"0-5, 7, 12-23, 25, 30-41, 43, 48-53"',\
            str(binsize),str(minPI), str(maxPI))
    print(lc_text)
    os.system(lc_text)
    return

def lelcgen(evtfile, outpath, binsize=1, minPI=0, maxPI=1035):
    lc_text = 'lelcgen evtfile=%s outfile=%s userdetid=%s eventtype=1 binsize=%s minPI=%s maxPI=%s clobber=yes'\
            %(evtfile[0], outpath+'le_lc_smfov',\
            '"0,2-4,6-10,12,14,20,22-26,28-30,32,34-36,38-42,44,46,52,54-58,60-62,64,66-68,70-74,76,78,84,86-90,92-94"',\
            str(binsize),str(minPI), str(maxPI))
    print lc_text
    os.system(lc_text)
    return

def mk_lcdir(outpath):
    if not os.path.isdir(outpath):
        os.system('mkdir -p '+outpath)

def genlist(outpath, instru='HE'):
    if instru == 'HE':
        listfile = [os.path.join(product_path,instru,'lightcurve','he_lc_g0_0-17.lc')]
    elif instru == 'ME':
        listfile = [os.path.join(product_path,instru,'lightcurve','me_lc_smfov_g0_0-53.lc')]
    else:
        listfile = [os.path.join(product_path,instru,'lightcurve','le_lc_smfov_g0_0-94.lc')]
    outname = os.path.join(outpath,'lightcurve.lst')
    with open(outname,'w') as fout:
        for item in listfile:
            wrt_str = item + '\n'
            fout.write(wrt_str)
    return outname

def hebkgmap(blindfile, ehkfile, gtifile, deadfile, listfile, outpath, minPI=20, maxPI=200):
    lcbkgmap_text = 'hebkgmap lc %s %s %s %s %s %s %s %s'\
            %(blindfile[0], ehkfile[0], gtifile[0], deadfile[0], listfile, str(minPI), str(maxPI),\
            os.path.join(outpath,'lc_bkgmap'))
    print lcbkgmap_text
    os.system(lcbkgmap_text)
    return

def mebkgmap(blindfile, ehkfile, gtifile, deadfile, tempfile, listfile, outpath, minPI=0, maxPI=1024):
    lcbkgmap_text = 'mebkgmap lc %s %s %s %s %s %s %s %s %s'\
            %(blindfile[0], ehkfile[0], gtifile[0], deadfile[0], tempfile[0], listfile, str(minPI), str(maxPI),\
            os.path.join(outpath,'lc_bkgmap'))
    print lcbkgmap_text
    os.system(lcbkgmap_text)
    return
    
def lebkgmap(blindfile, gtifile, listfile, outpath, minPI=0, maxPI=1535):
    lcbkgmap_text = 'lebkgmap lc %s %s %s %s %s %s'\
            %(blindfile[0], gtifile[0], listfile, str(minPI), str(maxPI),\
            os.path.join(outpath,'lc_bkgmap'))
    print lcbkgmap_text
    os.system(lcbkgmap_text)
    return

if __name__ == '__main__':
    outpath = os.path.join(product_path,'HE','lightcurve/')
    mk_lcdir(outpath)
    screenfile, blindfile, ehkfile, gtifile, _, deadfile = readfile(instru='HE')
    helcgen(screenfile, outpath, deadfile)
    listfile = genlist(outpath, instru='HE')
    hebkgmap(blindfile, ehkfile, gtifile, deadfile, listfile, outpath, minPI=20, maxPI=200)

    outpath = os.path.join(product_path,'ME','lightcurve/')
    mk_lcdir(outpath)
    screenfile, blindfile, ehkfile, gtifile, tempfile, deadfile = readfile(instru='ME')
    melcgen(screenfile, outpath, deadfile)
    mebkgmap(blindfile, ehkfile, gtifile, deadfile, tempfile, listfile, outpath, minPI=0, maxPI=1024)

    outpath = os.path.join(product_path,'LE','lightcurve/')
    mk_lcdir(outpath)
    screenfile, blindfile, ehkfile, gtifile, tempfile, deadfile = readfile(instru='LE')
    lelcgen(screenfile, outpath)
    lebkgmap(blindfile, gtifile, listfile, outpath, minPI=0, maxPI=1535)



