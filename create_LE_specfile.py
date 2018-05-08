#!/usr/bin/env python
from __future__ import division
from astropy.io import fits
import numpy as np

Channel = np.linspace(0,1535,1536)
counts  = np.zeros(len(Channel))
quality = np.zeros(len(Channel))
grouping = np.zeros(len(Channel))
exposure = 699.6667
outfile = '/home/tuoyl/Test/le_res_spec_test/new2.fits'

# Table
c1 = fits.Column(name='Channel', array=Channel,format='1J')
c2 = fits.Column(name='COUNTS' , array=counts ,format='1J')
c3 = fits.Column(name='QUALITY', array=quality,format='1I')
c4 = fits.Column(name='GROUPING', array=grouping,format='1I')
tb = fits.BinTableHDU.from_columns([c1,c2,c3,c4])
# Prmary Header
header = fits.Header()
primary_hdr = fits.Header()
primary_hdr['comments'] = 'FITS(Flexible Image Transport System)'
primary_hdu = fits.PrimaryHDU(header=primary_hdr)

hdul = fits.HDUList([primary_hdu, tb])
hdul.writeto(outfile,overwrite=True)

# write Keywords
TRUE = np.bool(True)
FALSE = np.bool(False)
hdulist = fits.open(outfile)
hdulist[1].header['EXTNAME'] = 'SPECTRUM'
hdulist[1].header['DETCHANS'] = 1536
hdulist[1].header['BACKFILE'] = 'NONE'
hdulist[1].header['BACKSCAL'] = 1
hdulist[1].header['CORRFILE'] = 'NONE'
hdulist[1].header['CORRSCAL'] = 0
hdulist[1].header['RESPFILE'] = 'NONE'
hdulist[1].header['ANCRFILE'] = 'NONE'
hdulist[1].header['FILETER']  = 'NONE'
hdulist[1].header['PHAVERSN'] = '1992a'
hdulist[1].header['STATERR']  = FALSE
hdulist[1].header['SYSERR']   = FALSE
hdulist[1].header['POISSERR'] = TRUE
hdulist[1].header['GROUPING'] = 1
hdulist[1].header['QUALITY']  = 1
hdulist[1].header['AREASCAL'] = 1
hdulist[1].header['EXPOSURE'] = exposure
hdulist[1].header['LIVETIME'] = 1
hdulist[1].header['DEADC']    = 0
hdulist[1].header['DETID']    = 0
hdulist[1].header['CHANTYPE'] = 'PI'
hdulist[1].header['TLMIN2']   = 0
hdulist[1].header['TLMAX2']   = 1535
hdulist[1].header['TELESCOP'] = 'HXMT'
hdulist[1].header['INSTRUME'] = 'LE'


hdulist.writeto(outfile,overwrite=True)
