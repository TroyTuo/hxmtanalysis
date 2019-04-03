#!/usr/bin/env python
from __future__ import division
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys


class InputRuleError(Exception):
    """raised when the input text file contained different lines of content"""
    pass

class IOError(Exception):
    """raised when the number of input arguments is not 6,
    which containing spectrum list, background list, response list,
    and their output file name"""
    pass



def usage():
    print("USAGE:")
    print("      ./addspec.py spec.txt bkg.txt rsp.txt outspec.pha outbkg.pha outrsp.fits\n")
    print("spec.txt: The input is an ASCII file containing a list of PHA datasets, one line for one file")
    print("bkg.txt:  The input is an ASCII file containing a list of background datasets, one line for one file")
    print("rsp.txt:  The input is an ASCII file containing a list of response datasets, one line for one file\n")
    print("EXAMPLE:")
    print("      if you do not need to merge spectra, background or response file, set the corresponding file to 'NONE' or 'none\n")
    print("      ./addspec.py spec.txt NONE NONE outspec.pha NONE NONE")
    sys.exit()


def read_input():
    if sys.argv[1] == "-h" or sys.argv[1] == "-H":
        usage()
    if len(sys.argv) is not 7:
        raise IOError("Input error, try -h to see the usage")
    if sys.argv[1] == "none" or sys.argv[1] == "NONE":
        specfile = False
    else:
        specfile = sys.argv[1]

    if sys.argv[2] == "none" or sys.argv[2] == "NONE":
        backfile = False
    else:
        backfile = sys.argv[2]

    if sys.argv[3] == "none" or sys.argv[3] == "NONE":
        rspfile = False
    else:
        rspfile = sys.argv[3]

    if sys.argv[4] == "none" or sys.argv[4] == "NONE":
        outspec = False
    else:
        outspec = sys.argv[4]

    if sys.argv[5] == "none" or sys.argv[5] == "NONE":
        outback = False
    else:
        outback = sys.argv[5]

    if sys.argv[6] == "none" or sys.argv[6] == "NONE":
        outrsp = False
    else:
        outrsp = sys.argv[6]

    return specfile, backfile, rspfile, outspec, outback, outrsp

def read_files(specfile, backfile, rspfile):
    if specfile:
        speclist = open(specfile)
        spectra = speclist.readlines()
        spectra = [x[0:-1] for x in spectra]
    else:
        spectra = False

    if backfile:
        backlist = open(backfile)
        background = backlist.readlines()
        background = [x[0:-1] for x in background]
    else:
        background = False

    if rspfile:
        rsplist = open(rspfile)
        response = rsplist.readlines()
        response = [x[0:-1] for x in response]
    else:
        response = False

    return spectra, background, response

def get_fits_info(filename):
    hdulist = fits.open(filename)
    exposure = np.float(hdulist[1].header["EXPOSURE"])
    instrument = hdulist[1].header["INSTRUME"]
    telescop   = hdulist[1].header["TELESCOP"]
    backfile   = hdulist[1].header['BACKFILE']
    chanum     = np.int(hdulist[1].header['DETCHANS'])
    respfile   = hdulist[1].header['RESPFILE']

    return exposure, instrument, telescop, backfile, chanum, respfile

def merge_spec(spectra, outspec):

    if not spectra or not outspec:
        return

    exposure_tmp, instrument_tmp, telescop_tmp, backfile_tmp, chanum_tmp, respfile_tmp = get_fits_info(spectra[0])

    expsure_new = 0.0
    instrument_new = instrument_tmp
    counts_new = np.zeros(chanum_tmp)


    for spectrum in spectra:

        exposure, instrument, telescop, backfile, chanum, respfile = get_fits_info(spectrum)

        if instrument == instrument_new:
            instrume_new = instrument
            exposure_new = expsure_new + exposure

            hdulist = fits.open(spectrum)
            counts  = hdulist[1].data.field("COUNTS")
            counts_new = counts_new + counts
        else:
            print("ERROR: the spectra were not detected by same INSTRUMENT!")
            sys.exit()

    create_specfile(chanum, instrument_new, exposure_new, counts_new, outspec, backfile='NONE', respfile='NONE', telescop=telescop)
    return 

def merge_bkg(backgrounds, outspec, spectra='', mission='HXMT'):

    if not backgrounds or not outspec:
        return

    if mission is not "HXMT":
        merge_spec(backgrounds, outspec)
        return

    exposure_spectmp, instrument_spectmp, telescop_spectmp, backfile_spectmp, chanum_spectmp, respfile_spectmp = get_fits_info(spectra[0])
    exposure_bkg, instrument_bkg, telescop_bkg, backfile_bkg, chanum_bkg, respfile_bkg = get_fits_info(backgrounds[0])

    expsure_new = 0.0
    instrument_new = instrument_bkg
    counts_new = np.zeros(chanum_bkg)


    for background, spectrum in zip(backgrounds, spectra):

        exposure_spec, instrument_spec, telescop_spec, backfile_spec, chanum_spec, respfile_spec = get_fits_info(spectrum)
        exposure, instrument, telescop, backfile, chanum, respfile = get_fits_info(background)

        if instrument == instrument_new:
            instrume_new = instrument
            #use spectrum exposure time
            exposure_new = expsure_new + exposure_spec

            hdulist = fits.open(background)
            counts  = hdulist[1].data.field("COUNTS")
            # weight by exposure
            counts  = (counts * exposure_spec)/exposure
            counts_new = counts_new + counts
        else:
            print("ERROR: the spectra were not detected by same INSTRUMENT!")
            sys.exit()

    create_specfile(chanum, instrument_new, exposure_new, counts_new, outspec, backfile='NONE', respfile='NONE', telescop=telescop)
    return 


def merge_rsp(rspfiles, outrsp, specfiles, bkgfiles):

    if not rspfiles or not outrsp:
        return
    if not specfiles or not bkgfiles:
        print("WARNING: NO spectra file and background were found, the weight for each rsp is the same")

    #NOTE:Only the Matrix were calculated, the EBoundary were copied from old resp
    hdulist = fits.open(rspfiles[0])
    ENERG_LO = hdulist[1].data.field("ENERG_LO")
    ENERG_HI = hdulist[1].data.field("ENERG_HI")
    N_GRP    = hdulist[1].data.field("N_GRP")
    F_CHAN   = hdulist[1].data.field("F_CHAN")
    N_CHAN   = hdulist[1].data.field("N_CHAN")

    CHANNEL  = hdulist[2].data.field("CHANNEL")
    E_MIN    = hdulist[2].data.field("E_MIN")
    E_MAX    = hdulist[2].data.field("E_MAX")

    TELESCOP = hdulist[1].header['TELESCOP']
    INSTRUME  = hdulist[1].header['INSTRUME']
    DETCHANS = hdulist[1].header['DETCHANS']

    matrix_tmp = []
    exposure_tmp = 0.0
    for resp, spectrum, background in zip(rspfiles, specfiles, bkgfiles):
        exposure_spec, _, _, _, _, _ = get_fits_info(spectrum)
        exposure_bkg , _, _, _, _, _ = get_fits_info(background)
        hdulist = fits.open(resp)
        matrix  = hdulist[1].data.field("MATRIX")
        matrix_tmp.append(matrix * exposure_spec)
        exposure_tmp = exposure_tmp + exposure_spec
        matrix_format = hdulist[1].header['TFORM6']
    matrix_i = matrix_tmp[0]
    for i in xrange(len(matrix_tmp)):
        if i == 0 :continue
        matrix_i = matrix_i + matrix_tmp[i]
    matrix = matrix_i / exposure_tmp
    create_rspfile(ENERG_LO, ENERG_HI, N_GRP, F_CHAN, N_CHAN, matrix, 
            CHANNEL, E_MIN, E_MAX, matrix_format, outrsp,
            telescop=TELESCOP, instrume=INSTRUME, detchans=DETCHANS)



def create_rspfile(energ_lo, energ_hi, n_grp, f_chan, n_chan, matrix, 
    channel, e_min, e_max, matrix_format, outfile,
    telescop='HXMT', instrume='HE', detchans=256):

    # Table
    c1 = fits.Column(name='ENERG_LO', array=energ_lo,format='1E')
    c2 = fits.Column(name='ENERG_HI', array=energ_hi ,format='1E')
    c3 = fits.Column(name='N_GRP',  array=n_grp,format='1I')
    c4 = fits.Column(name='F_CHAN', array=f_chan,format='1I')
    c5 = fits.Column(name='N_CHAN', array=n_chan,format='1I')
    c6 = fits.Column(name='MATRIX', array=matrix,format=matrix_format)

    cc1 = fits.Column(name='CHANNEL', array=channel,format='1I')
    cc2 = fits.Column(name='E_MIN',   array=e_min,format='1E')
    cc3 = fits.Column(name='E_MAX',   array=e_max,format='1E')

    tb1 = fits.BinTableHDU.from_columns([c1,c2,c3,c4,c5,c6])
    tb2 = fits.BinTableHDU.from_columns([cc1,cc2,cc3])
    # Prmary Header
    header = fits.Header()
    primary_hdr = fits.Header()
    primary_hdr['comments'] = 'FITS(Flexible Image Transport System)'
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    
    hdul = fits.HDUList([primary_hdu, tb1, tb2])
    hdul.writeto(outfile,overwrite=True)
    
    # write Keywords
    TRUE = np.bool(True)
    FALSE = np.bool(False)

    hdulist = fits.open(outfile)
    hdulist[1].header['TUNIT1']  = 'keV'
    hdulist[1].header['TUNIT2']  = 'keV'
    hdulist[1].header['EXTNAME'] = 'MATRIX'
    hdulist[1].header['TELESCOP']= telescop
    hdulist[1].header['INSTRUME']= instrume
    hdulist[1].header['DETNAM']  = instrume
    hdulist[1].header['DETCHANS']= detchans
    hdulist[1].header['CHANTYPE']= 'PI'
    hdulist[1].header['HDUVERS'] = '1.3.0'
    hdulist[1].header['CCLS0001']= 'BCF'
    hdulist[1].header['HDUCLASS']= 'OGIP'
    hdulist[1].header['HDUCLAS1']= 'RESPONSE'
    hdulist[1].header['HDUCLAS2']= 'RSP_MATRIX'
    hdulist[1].header['CDTP0001']= 'DATA'
    hdulist[1].header['TLMIN4']  = '0'

    hdulist[2].header['TUNIT2']  = 'keV'
    hdulist[2].header['TUNIT3']  = 'keV'
    hdulist[2].header['EXTNAME'] = 'EBOUNDS'
    hdulist[2].header['TELESCOP']= telescop
    hdulist[2].header['INSTRUME']= instrume
    hdulist[2].header['DETNAM']  = instrume
    hdulist[2].header['DETCHANS']= detchans
    hdulist[2].header['CHANTYPE']= 'PI'
    hdulist[2].header['HDUVERS'] = '1.3.0'
    hdulist[2].header['CCLS0001']= 'BCF'
    hdulist[2].header['HDUCLASS']= 'OGIP'
    hdulist[2].header['HDUCLAS1']= 'RESPONSE'
    hdulist[2].header['HDUCLAS2']= 'EBOUNDS'
    hdulist[2].header['CDTP0001']= 'DATA'
    hdulist[2].header['TLMIN4']  = '0'


    hdulist.writeto(outfile,overwrite=True)

    return


def create_specfile(chanum, instrume, exposure, counts, outfile, backfile='NONE', respfile='NONE', telescop='HXMT'):
    Channel = np.linspace(0,chanum-1,chanum)
    quality = np.zeros(len(Channel))
    COUNTS = counts
    grouping = np.zeros(len(Channel))
    
    # Table
    c1 = fits.Column(name='Channel', array=Channel,format='1J')
    c2 = fits.Column(name='COUNTS' , array=COUNTS ,format='1J')
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
    hdulist[1].header['DETCHANS'] = chanum
    hdulist[1].header['BACKFILE'] = backfile
    hdulist[1].header['BACKSCAL'] = 1
    hdulist[1].header['CORRFILE'] = 'NONE'
    hdulist[1].header['CORRSCAL'] = 0
    hdulist[1].header['RESPFILE'] = respfile
    hdulist[1].header['ANCRFILE'] = 'NONE'
    hdulist[1].header['FILETER']  = 'NONE'
    hdulist[1].header['PHAVERSN'] = '1992a'
    hdulist[1].header['STATERR']  = FALSE
    hdulist[1].header['SYSERR']   = FALSE
    hdulist[1].header['POISSERR'] = TRUE
    #NOTE:the error of spetrum is the poisson error of counts, no error assigned.
    hdulist[1].header['GROUPING'] = 1
    hdulist[1].header['QUALITY']  = 1
    hdulist[1].header['AREASCAL'] = 1
    hdulist[1].header['EXPOSURE'] = exposure
    hdulist[1].header['LIVETIME'] = 1
    hdulist[1].header['DEADC']    = 0
    hdulist[1].header['DETID']    = 0
    hdulist[1].header['CHANTYPE'] = 'PI'
    hdulist[1].header['TLMIN2']   = 0
    hdulist[1].header['TLMAX2']   = chanum-1
    hdulist[1].header['TELESCOP'] = telescop
    hdulist[1].header['INSTRUME'] = instrume
    
    hdulist.writeto(outfile,overwrite=True)

if __name__ == "__main__":
    specfile, backfile, rspfile, outspec, outback, outrsp = read_input()
    spectra, background, response = read_files(specfile, backfile, rspfile)
    merge_spec(spectra, outspec)
    merge_bkg(background, outback, spectra, mission='HXMT')
    merge_rsp(response, outrsp, spectra, background)
    print("DONE")
