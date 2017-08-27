from function import *

raw_filename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-Evt_FFFFFF_V1_1RP.FITS'
raw_filename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/Event4piV2_select.fits'

pi_filename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/hepicalV2.fits'
pi_filename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/hepicalV2_selected.fits'

raw_data,det_id,channel,pulse_width,event_type = readdata(raw_filename)
pi_data,pi_det_id,pi_channel,pi_pulse_width,pi_event_type = readdata(pi_filename)




lc_time1,lc1 = genlc(raw_data,binsize = 0.01)
#plt.figure()
#plt.step(lc_time1,lc1,color='blue')
#lc_time2,lc2 = genlc(pi_data,binsize = 0.01)
#plt.step(lc_time2,lc2,color='green')
#plt.title('GRS 1915')
#plt.ylabel('rate(cnt/s)')
#plt.xlabel('time(s)')
#
#how much event_type=0 and event_type = 1 events are expelled 
evt0_index = np.where(event_type == 0)
evt0_evt = raw_data[evt0_index]
evt0_lc_time,evt0_lc = genlc(evt0_evt,binsize=0.01)

evt1_index = np.where(event_type == 1)
evt1_evt = raw_data[evt1_index]
evt1_lc_time,evt1_lc = genlc(evt1_evt,binsize=0.01)

pi_evt0_index = np.where(pi_event_type == 0)
pi_evt0_evt = pi_data[pi_evt0_index]
pi_evt0_lc_time,pi_evt0_lc = genlc(pi_evt0_evt,binsize=0.01)

pi_evt1_index = np.where(pi_event_type == 1)
pi_evt1_evt = pi_data[pi_evt1_index]
pi_evt1_lc_time,pi_evt1_lc = genlc(pi_evt1_evt,binsize=0.01)

plt.figure('evt_type')
plt.subplot(4,1,1)
plt.title('evt_type0')
plt.step(evt0_lc_time,evt0_lc,color='blue')
plt.step(pi_evt0_lc_time,pi_evt0_lc,color='green')

plt.subplot(4,1,2)
plt.title('ratio')
plt.step(evt0_lc_time,evt0_lc/pi_evt0_lc)

plt.subplot(4,1,3)
plt.title('evt_type1')
plt.step(evt1_lc_time,evt1_lc,color='blue')
plt.step(pi_evt1_lc_time,pi_evt1_lc,color='green')

plt.subplot(4,1,4)
plt.title('ratio')
plt.step(evt1_lc_time,evt1_lc/pi_evt1_lc)

# time intervals analysis (for physical events of one detector)

## delta time of whole detected photon of rawdata
#plt.figure('all_data')
#ti_all = []
#tmp_raw_data = np.sort(raw_data)
#for i in xrange(len(raw_data)-1):
#    ti_all.append(tmp_raw_data[i+1]-tmp_raw_data[i])
#all_ti_x,all_ti_y = genlc(ti_all,binsize=0.00001,rate=False)
#plt.semilogy(1e6*all_ti_x,all_ti_y)
#plt.ylabel('cnt')
#plt.xlabel('delta_time($\mu s$)')
#plt.title('all_data')
#plt.xlim([-100,1550])
#plt.ylim([0,1e6])
#
#
#det_num = 0 # select detector 0 
#event_type_num = 0
#
#    #select detector 
#det_id_index = np.where(det_id == det_num)
#event_type_index = np.where(event_type == event_type_num)
#
#    #select event type
#pi_det_id_index = np.where(pi_det_id == det_num)
#pi_event_type_index = np.where(pi_event_type == event_type_num)
#
#select_index = np.concatenate((det_id_index,event_type_index),axis=1)
#select_index = np.unique(select_index)
#print select_index
#pi_select_index = np.concatenate((pi_det_id_index,pi_event_type_index),axis=1)
#pi_select_index = np.unique(pi_select_index)
#print pi_select_index
#print len(raw_data),len(pi_data)
#
#select_index = det_id_index
#pi_select_index = pi_det_id_index
#print len(select_index),len(pi_select_index)
#
#select_evt = raw_data[select_index]
#pi_select_evt = pi_data[pi_select_index]
#select_evt.sort()
#pi_select_evt.sort()
#print len(select_evt),len(pi_select_evt)
#ti = []
#pi_ti = []
#for i in xrange(len(select_evt)-1):
#    ti.append(select_evt[i+1]-select_evt[i])
#for i in xrange(len(pi_select_evt)-1):
#    pi_ti.append(pi_select_evt[i+1]-pi_select_evt[i])
#
#print len(ti),len(pi_ti)
#
#plt.figure('delta_time')
#plt.subplot(2,1,1)
#plt.title('delta_time')
#ti_x,ti_y = genlc(ti,binsize=0.00001,rate=False)
#plt.semilogy(1e6*ti_x,ti_y)
#plt.xlim([-1000,15000])
#plt.ylim([1e0,1e4])
#print "min time_interval for rawdata: ",float(min(ti))
#plt.ylabel('cnt')
#plt.subplot(2,1,2)
#pi_ti_x,pi_ti_y = genlc(pi_ti,binsize=0.00001,rate=False)
#plt.semilogy(1e6*pi_ti_x,pi_ti_y)
#plt.ylabel('cnt')
#plt.subplot(2,1,2)
#plt.xlabel('delta_time($\mu s$)')
#plt.xlim([-1000,15000])
#plt.ylim([1e0,1e4])
#print "min time_interval for pidata: ",float(min(pi_ti))
#

# select one PUCD, AKA det 0-5
det_id0_index = np.where(det_id == 0)
det_id1_index = np.where(det_id == 1)
det_id2_index = np.where(det_id == 2)
det_id3_index = np.where(det_id == 3)
det_id4_index = np.where(det_id == 4)
det_id5_index = np.where(det_id == 5)
det_id_index = np.concatenate((det_id0_index,det_id1_index,det_id2_index,det_id3_index,det_id4_index,det_id5_index),axis=1)
det_id_index = np.unique(det_id_index)

pi_det_id0_index = np.where(pi_det_id == 0)
pi_det_id1_index = np.where(pi_det_id == 1)
pi_det_id2_index = np.where(pi_det_id == 2)
pi_det_id3_index = np.where(pi_det_id == 3)
pi_det_id4_index = np.where(pi_det_id == 4)
pi_det_id5_index = np.where(pi_det_id == 5)
pi_det_id_index = np.concatenate((pi_det_id0_index,pi_det_id1_index,pi_det_id2_index,pi_det_id3_index,pi_det_id4_index,pi_det_id5_index),axis=1)
pi_det_id_index = np.unique(pi_det_id_index)

evt_05 = raw_data[det_id_index]
evt_05.sort()

pi_evt_05 = pi_data[pi_det_id_index]
pi_evt_05.sort()
print 'lllllllllllll',evt_05,pi_evt_05
print len(evt_05),len(pi_evt_05)
ti05 = []
pi_ti05 = []
for i in xrange(len(evt_05)-1):
    ti05.append(evt_05[i+1]-evt_05[i])
for i in xrange(len(pi_evt_05)-1):
    pi_ti05.append(pi_evt_05[i+1]-pi_evt_05[i])
plt.figure('det_0-5')
plt.title('det_0-5')
plt.subplot(2,1,1)
lc_x,lc_y = genlc(ti05,binsize=0.00001,rate=False)
plt.semilogy(1e6*lc_x,lc_y)
plt.subplot(2,1,2)
pi_lc_x,pi_lc_y = genlc(pi_ti05,binsize=0.00001,rate=False)
plt.semilogy(1e6*pi_lc_x,pi_lc_y)




#select those spikes data

 


#lc_x,lc_y = genlc(raw_data,binsize=0.001)
#plt.figure()
#plt.step(lc_x,lc_y)




plt.show()
