#! /bin/sh
#################
# Tuoyl
from function import *

datafilename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-Evt_FFFFFF_V1_1RP.FITS'\

hvfilename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-HV_FFFFFF_V1_1RP.FITS'

pmfilename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-PM_FFFFFF_V1_1RP.FITS'

ehkfilename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/AUX/EHK_man.fits'

saafilename = '/hxmt/home/hxmtsoft/hxmtehkgen/hxmtehkgen/SAA/SAA.fits'

raw_data,det_id,channel,pulse_width,event_type = readdata(datafilename,sat='hxmt')

# SAA analysis
hdulist_saa = pf.open(saafilename)
tb = hdulist_saa[1].data
saa_start = tb.field(0)
saa_stop = tb.field(1)



# HV file analysis

    # light curves
bin = 1
lc_time,lc = genlc(raw_data,binsize=bin)
lc_rate = lc/bin
plt.subplot(7,1,1)
plt.ylabel('rate(cnt/s)')
plt.plot(lc_time,lc_rate)
lc_mincut = lc_time[min(np.where( lc_rate == 0)[0])]
lc_maxcut = lc_time[max(np.where( lc_rate == 0)[0])]
print "lc_min_cut",lc_mincut
print "lc_man_cut",lc_maxcut

saa_mincut = []
saa_maxcut = []
for i in xrange(len(saa_start)):
    if saa_start[i] < max(lc_time) and saa_start[i] > min(lc_time):
        saa_mincut.append(saa_start[i])
        saa_maxcut.append(saa_stop[i])

plt.axvline(x = lc_mincut,color='black')
plt.axvline(x = lc_maxcut,color='black')
plt.axvline(x = saa_mincut[0],color='red')
plt.axvline(x = saa_maxcut[0],color='red')




    # HV phowise detectors
hdulist = pf.open(hvfilename)
tb = hdulist[1].data
time_hv_ph = tb.field(0)
hv_ph0 = tb.field(1)
#hv_ph1 = tb.field(2)
#hv_ph2 = tb.field(3)
#hv_ph3 = tb.field(4)
#hv_ph4 = tb.field(5)
#hv_ph5 = tb.field(6)
#hv_ph6 = tb.field(7)
#hv_ph7 = tb.field(8)
#hv_ph8 = tb.field(9)
#hv_ph9 = tb.field(10)
#hv_ph10 = tb.field(11)
#hv_ph11 = tb.field(12)
#hv_ph12 = tb.field(13)
#hv_ph13= tb.field(14)
#hv_ph14 = tb.field(15)
#hv_ph15 = tb.field(16)
#hv_ph16 = tb.field(17)
#hv_ph17 = tb.field(18)

plt.subplot(7,1,2)
plt.title('hv_ph')
plt.xlabel('time')
plt.ylabel('Volt')

for i in xrange(len(time_hv_ph)-1):
    if abs(time_hv_ph[i] - time_hv_ph[i+1]) > 1:
        hv_ph_mincut = time_hv_ph[i]
        hv_ph_maxcut = time_hv_ph[i+1]
print "hv_ph_mincut",hv_ph_mincut
print "hv_ph_maxcut",hv_ph_maxcut
plt.axvline(x=hv_ph_mincut,color='black')
plt.axvline(x=hv_ph_maxcut,color='black')

plt.axvline(x = saa_mincut[0],color='red')
plt.axvline(x = saa_maxcut[0],color='red')

plt.plot(time_hv_ph,hv_ph0)
#plt.plot(time_hv_ph,hv_ph1)
#plt.plot(time_hv_ph,hv_ph2)
#plt.plot(time_hv_ph,hv_ph3)
#plt.plot(time_hv_ph,hv_ph4)
#plt.plot(time_hv_ph,hv_ph5)
#plt.plot(time_hv_ph,hv_ph6)
#plt.plot(time_hv_ph,hv_ph7)
#plt.plot(time_hv_ph,hv_ph8)
#plt.plot(time_hv_ph,hv_ph9)
#plt.plot(time_hv_ph,hv_ph10)
#plt.plot(time_hv_ph,hv_ph11)
#plt.plot(time_hv_ph,hv_ph12)
#plt.plot(time_hv_ph,hv_ph13)
#plt.plot(time_hv_ph,hv_ph14)
#plt.plot(time_hv_ph,hv_ph15)
#plt.plot(time_hv_ph,hv_ph16)
#plt.plot(time_hv_ph,hv_ph17)
    
    # HV occ detector
tb1 = hdulist[2].data
time_hv_ooc = tb1.field(0)
hv_ooc0 = tb1.field(1)
hv_ooc1 = tb1.field(2)
hv_ooc2 = tb1.field(3)
plt.subplot(7,1,3)
plt.title('hv_occ')
plt.ylabel('Volt')

for i in xrange(len(time_hv_ooc)-1):
    if abs(time_hv_ooc[i] - time_hv_ooc[i+1]) > 1:
        hv_ooc_mincut = time_hv_ooc[i]
        hv_ooc_maxcut = time_hv_ooc[i+1]
print "hv_ooc_mincut",hv_ooc_mincut
print "hv_ooc_maxcut",hv_ooc_maxcut
plt.axvline(x=hv_ooc_mincut,color='black')
plt.axvline(x=hv_ooc_maxcut,color='black')
plt.axvline(x = saa_mincut[0],color='red')
plt.axvline(x = saa_maxcut[0],color='red')
plt.plot(time_hv_ooc,hv_ooc0)
#plt.plot(time_hv_ooc,hv_ooc1)
#plt.plot(time_hv_ooc,hv_ooc2)

    # HV PM detector
tb2 = hdulist[3].data
time_hv_pm = tb1.field(0)
hv_pm0 = tb1.field(1)
hv_pm1 = tb1.field(2)
hv_pm2 = tb1.field(3)
plt.subplot(7,1,4)
plt.title('hv_pm')
plt.ylabel('Volt')
for i in xrange(len(time_hv_pm)-1):
    if abs(time_hv_pm[i] - time_hv_pm[i+1]) > 1:
        hv_pm_mincut = time_hv_pm[i]
        hv_pm_maxcut = time_hv_pm[i+1]
print "hv_pm_mincut",hv_pm_mincut
print "hv_pm_maxcut",hv_pm_maxcut
plt.axvline(x=hv_pm_mincut)
plt.axvline(x=hv_pm_maxcut)
plt.plot(time_hv_pm,hv_pm0)
#plt.plot(time_hv_pm,hv_pm1)
#plt.plot(time_hv_pm,hv_pm2)

hdulist.close()


# PM counts analysis
hdulist_pm = pf.open(pmfilename)
tb = hdulist_pm[1].data
time_pm = tb.field(0)
cnt_pm_0 = tb.field(1)
cnt_pm_1 = tb.field(2)
cnt_pm_2 = tb.field(3)
hdulist_pm.close()

plt.subplot(7,1,5)
plt.plot(time_pm,cnt_pm_0)
#plt.plot(time_pm,cnt_pm_1)
#plt.plot(time_pm,cnt_pm_2)
plt.title('cnt_pm')
plt.ylabel('cnt/s')

plt.axvline(x = saa_mincut[0],color='red')
plt.axvline(x = saa_maxcut[0],color='red')


# EHK analysis
hdulist_ehk = pf.open(ehkfilename)
tb = hdulist_ehk[1].data
time_ehk = tb.field(0)
saa_flag = tb.field(22)
hdulist_ehk.close()

plt.subplot(7,1,6)
plt.plot(time_ehk,saa_flag)
plt.title('saa_flag')

# SAA analysis
hdulist_saa = pf.open(saafilename)
tb = hdulist_saa[1].data
saa_start = tb.field(0)
saa_stop = tb.field(1)


plt.show()


































