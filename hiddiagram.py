from function import *
def readscreen(filename):
    hdulist = pf.open(filename)
    tb = hdulist[1].data
    time = tb.field(0)
    pi = tb.field(1)
    det_id = tb.field(2)
    hdulist.close()
    return time,pi,det_id

ratio = []
intensity = []
evt_all = []
time = []

#read data
with open('/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/lightcurve_all/screen.lst','r') as f:
    read_filename = f.read()
    line = read_filename.split('\n')
    
for filename in line:
    if filename == '':continue
    raw_data,channel,det_id = readscreen(filename)
    if (175010114 >= min(raw_data) and 175010114 <= max(raw_data))or (175011504 <= max(raw_data) and 175011504 >= min(raw_data)):continue
    
    evt_all.append(raw_data)
    # parameters for HID cut
    channel_cut1 = 150 
    channel_cut2 =150
    time_set = 20 # how much time duration (s) for one dot on HID
    N =int((max(raw_data)-min(raw_data))/time_set)
    evt_dur = len(channel)/N #events
    # generate time cut
    evt_index = []
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

if evt_all==[]:print "NO DATA"
plt.figure('HID diagram')
plt.title('HID diagram')
plt.xlabel('')
plt.loglog(intensity,ratio,linestyle='-', marker='o')

plt.figure('HID evolution')
plt.title('HID evolution')
evt_all = np.concatenate(evt_all)
lc_x, lc_y = genlc(evt_all,binsize=10,fig=False,rate=True)
plt.subplot(2,1,1)
plt.title('light curve')
plt.axvline(x=175010114,color='r')
plt.axvline(x=175011504,color='r')
plt.plot(lc_x,lc_y)
plt.subplot(2,1,2)
plt.title('hardness evolution')
plt.axvline(x=175010114,color='r')
plt.axvline(x=175011504,color='r')
plt.plot(time,ratio,linestyle='-',marker='+')

plt.show()
