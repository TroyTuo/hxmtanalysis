#!/usr/bin/env python
#coding: utf-8

'''
use bkgmap python script to generate background for spec and lc
'''

import argparse
import os
import glob
from natsort import natsorted, ns
from astropy.io import fits


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
        specfile = glob.glob(os.path.join(product_path,instru,'spectra','he_spec_NaI*.pha*'))
        specfile = sorthespec(specfile)
        genlist(specfile, os.path.join(product_path,instru,'spectra','spec.lst'))
    elif instru == 'ME':
        deadfile = glob.glob(os.path.join(product_path,instru,'tmp','*dtime*'))
        specfile = glob.glob(os.path.join(product_path,instru,'spectra','*smfov*.pha*'))
        genlist(specfile, os.path.join(product_path,instru,'spectra','spec.lst'))
    else:
        deadfile = []
        specfile = glob.glob(os.path.join(product_path,instru,'spectra','*smfov*.pha*'))
        genlist(specfile, os.path.join(product_path,instru,'spectra','spec.lst'))
    listfile = [os.path.join(product_path,instru,'spectra','spec.lst')]

    return blindfile, ehkfile, gtifile, tempfile, deadfile, specfile, listfile

def sorthespec(specfiles):
    new_specfiles = natsorted(specfiles, alg=ns.IGNORECASE)
    cp_specfiles = new_specfiles[1:-1]
    cp_specfiles.append(new_specfiles[0])
    cp_specfiles.append(new_specfiles[-1])
    new_specfiles = cp_specfiles
    return new_specfiles

def genlist(listname, outfile):
    with open(outfile,'w') as fout:
        for item in listname:
            wrt_str = item + '\n'
            fout.write(wrt_str)

def write_back(specfilename, bkgfilename):
    hdulist = fits.open(specfilename)
    hdulist[1].header['BACKFILE'] = bkgfilename
    hdulist.writeto(specfilename, overwrite=True)
    
def hebkgmap():
    blindfile, ehkfile, gtifile, tempfile, deadfile, specfilelist, listfile = readfile(instru='HE')
    specbkgmap_text = 'hebkgmap spec '+blindfile[0] +' '+ \
            ehkfile[0] + ' ' + gtifile[0] + ' ' +  deadfile[0] + ' ' + listfile[0] + ' 0 255 '+\
            os.path.join(product_path,'HE','spectra','he_bkg_spec')
    print specbkgmap_text
    os.system(specbkgmap_text)
    #sort
    bkgfilelist = glob.glob(product_path+'/HE/spectra/he_bkg_spec*')
    bkgfilelist = natsorted(bkgfilelist, alg=ns.IGNORECASE)
    for specfile_i, bkgfile_i in zip(specfilelist, bkgfilelist):
        print specfile_i, bkgfile_i
        write_back(specfile_i, bkgfile_i)
    return

def mebkgmap():
    blindfile, ehkfile, gtifile, tempfile, deadfile, specfile, listfile = readfile(instru='ME')
    specbkgmap_text = 'mebkgmap spec '+blindfile[0] +' '+ \
            ehkfile[0] + ' ' + gtifile[0] + ' ' +  deadfile[0] + ' ' + tempfile[0] +' '+\
            listfile[0] + ' 0 1024 '+ os.path.join(product_path,'ME','spectra','me_bkg_spec')
    print specbkgmap_text
    os.system(specbkgmap_text)
    
def lebkgmap():
    blindfile, _, gtifile, _, _, _, listfile = readfile(instru='LE')
    specbkgmap_text = 'lebkgmap spec '+blindfile[0] +' '+ \
            gtifile[0] + ' ' + listfile[0] + ' 0 1535 '+ os.path.join(product_path,'LE','spectra','le_bkg_spec')
    print specbkgmap_text
    os.system(specbkgmap_text)

if __name__ == "__main__":
    hebkgmap()
#    mebkgmap()
#    lebkgmap()


