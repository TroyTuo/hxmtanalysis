from function import *

ratio = []
intensity = []

with open('/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/GRS1915_allEvt','r') as f:
    read_filename = f.read()
    line = read_filename.split('\n')
    
for filename in line:
    if filename == '':continue
    raw_data,det_id,channel,pulse_width,event_type = readdata(filename)
    
    # parameters for HID cut
    channel_cut1 = 150
    channel_cut2 =150
    time_set = 100 #s
    N =int((max(raw_data)-min(raw_data))/time_set)
    evt_dur = len(channel)/N #events
    # generate time cut
    evt_index = []
    time = []
    for i in xrange(N):
        evt_index.append(evt_dur*i)
        time.append(raw_data[evt_dur*i])
    print evt_index
    print time
    
    for j in xrange(len(evt_index)):
        if j == len(evt_index)-1:
            channel_cut = channel[evt_index[j]:-1]
            duration = max(raw_data[evt_index[j]:-1]) - min (raw_data[evt_index[j]:-1])
        else:
            channel_cut = channel[evt_index[j]:evt_index[j+1]]
            duration = max(raw_data[evt_index[j]:evt_index[j+1]]) - min(raw_data[evt_index[j]:evt_index[j+1]])
        print duration
        channel_x,channel_y = genspec(channel_cut)
    
        # cut spectrum into two parts
        cnt1 = np.sum(channel_y[0:channel_cut1])
        cnt2 = np.sum(channel_y[channel_cut1:-1])
        ratio.append(cnt1/float(cnt2))
        intensity.append(np.sum(channel_y[0:-1])/duration)
        
plt.figure('HID diagram')
plt.title('HID diagram')
plt.xlabel('')
plt.loglog(ratio,intensity,linestyle='-', marker='o')
#plt.plot(ratio1,ratio2,'ro')

#plt.plot(time,ratio1,color='r',label='channel0-'+str(channel_cut1)+'/'+str(channel_cut1)+'-'+str(channel_cut2))
#plt.plot(time,ratio2,color='blue',label='channel'+str(channel_cut1)+'-'+str(channel_cut2)+'/'+str(channel_cut2)+'-255')
#plt.legend(loc='lower left')
plt.show()
