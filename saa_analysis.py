#! /bin/sh
from function import *

datafilename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-Evt_FFFFFF_V1_1RP.FITS'\

hvfilename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-HV_FFFFFF_V1_1RP.FITS'

pmfilenmae = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-PM_FFFFFF_V1_1RP.FITS'

raw_data,det_id,channel,pulse_width,event_type = readdata(datafilename,sat='hxmt')

hdulist = pf.open(hvfilename)
tb = hdulist[1].data
time_hv_ph = tb.field(0)
print time_hv_ph
hv_ph0 = tb.field(1)
print hv_ph0