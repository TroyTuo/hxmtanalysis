from function import *

filename = '/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/P0101310/P0101310001/P010131000102-20170727-01-01/HE/HXMT_P010131000102_HE-Evt_FFFFFF_V1_1RP.FITS'
raw_data,det_id,channel,pulse_width,event_type = readdata(filename)

# parameters for HID cut
channel_cut1 = 70
channel_cut2 =150
N = 11
evt_dur = len(channel)/N #events
# generate time cut
evt_index = []
time = []
for i in xrange(N):
    evt_index.append(evt_dur*i)
    time.append(raw_data[evt_dur*i])
print evt_index
print time

ratio1 = []
ratio2 = []
for j in xrange(len(evt_index)):
    if j == len(evt_index)-1:
        channel_cut = channel[evt_index[j]:-1]
    else:
        channel_cut = channel[evt_index[j]:evt_index[j+1]]
    print time[j]
    channel_x,channel_y = genspec(channel_cut)
    
    
    
    
    
    cnt1 = np.sum(channel_y[0:channel_cut1])
    cnt2 = np.sum(channel_y[channel_cut1:channel_cut2])
    cnt3 = np.sum(channel_y[channel_cut2:-1])
#    plt.plot(channel_x,channel_y)
    print cnt1/float(cnt2),cnt2/float(cnt3)
    ratio1.append(cnt1/float(cnt2))
    ratio2.append(cnt2/float(cnt3))
plt.figure('HID diagram')
plt.title('HID diagram')
plt.xlabel('time')

plt.plot(time,ratio1,color='r',label='channel0-'+str(channel_cut1)+'/'+str(channel_cut1)+'-'+str(channel_cut2))
plt.plot(time,ratio2,color='blue',label='channel'+str(channel_cut1)+'-'+str(channel_cut2)+'/'+str(channel_cut2)+'-255')
plt.legend(loc='lower left')
plt.show()
