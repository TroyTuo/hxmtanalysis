#!/usr/bin/env python

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
ratio_error = []
intensity_error = []
evt_all = []
time = []

#read data
with open('/home/hxmt/tuoyl/HXMT_process/archived_data/HE/PointSrc/GSR_1915/lightcurve_all/screen.lst','r') as f:
    read_filename = f.read()
    line = read_filename.split('\n')
    
for filename in line:
    ratio = []
    intensity = []
    ratio_error = []
    intensity_error = []
    evt_all = []
    time = []
    if filename == '':continue
    raw_data,channel,det_id = readscreen(filename)
    #if min(raw_data) < 174980000 or max(raw_data > 175010000):continue # time cut
    evt_all.append(raw_data)
    # parameters for HID cut
    channel_cut1 = 150 
    channel_cut2 =150
    time_set = 20 # how much time duration (s) for one dot on HID
    N =int((max(raw_data)-min(raw_data))/time_set)
    evt_dur = len(channel)/N #events
    # generate time cut
    
    test_evt_index = []
    time_list = np.linspace(min(raw_data),max(raw_data),int((max(raw_data)-min(raw_data))/time_set))
    test_cnt, test_time = np.histogram(raw_data,bins=time_list)
    for test in test_time:
        test_evt_index.append(np.where(raw_data == test))
    
    
    evt_index = []
    for i in xrange(N):
        evt_index.append(evt_dur*i)
        time.append(raw_data[evt_dur*i])
    #print evt_index
    #print time
    for j in xrange(len(evt_index)):
        if j == len(evt_index)-1:
            channel_cut = channel[evt_index[j]:-1]
            duration = max(raw_data[evt_index[j]:-1]) - min (raw_data[evt_index[j]:-1])
        else:
            channel_cut = channel[evt_index[j]:evt_index[j+1]]
            duration = max(raw_data[evt_index[j]:evt_index[j+1]]) - min(raw_data[evt_index[j]:evt_index[j+1]])
        channel_x,channel_y = genspec(channel_cut)
    
        # cut spectrum into two parts
        cnt1 = np.sum(channel_y[0:channel_cut1])
        cnt2 = np.sum(channel_y[channel_cut1:-1])
        if np.sum(channel_y/duration) < 200: continue
        ratio_tmp = cnt1/float(cnt2)
        intensity_tmp = np.sum(channel_y/duration)
        ratio.append(ratio_tmp)
        intensity.append(intensity_tmp)
        #calculate error
        cnt1_err = np.sqrt(cnt1)
        cnt2_err = np.sqrt(cnt2)
        intensity_error.append(np.sqrt(intensity_tmp))
        ratio_error.append(ratio_tmp * (np.sqrt(np.power(cnt1_err,2)/np.power(cnt1,2) + np.power(cnt2_err,2)/np.power(cnt2,2))))

    evt_all = np.concatenate(evt_all)
    plt.figure('HID diagram')
    plt.subplot(2,1,1)
    lc_x, lc_y = genlc(evt_all,binsize=10,fig=False,rate=True)
    plt.title('light curve')
    plt.ylabel('cnt/s')
    plt.step(lc_x,lc_y)
    
    #plt.subplot(3,1,2)
    #plt.plot(time,ratio,marker='+')
    #plt.ylabel('hardness([0-150]/[150-255]')
    
    plt.subplot(2,1,2)        
    plt.title('HID diagram')
    plt.xlabel('intensity(cnt/s)')
    plt.ylabel('hardness([0-150]/[150-255]')
    plt.errorbar(intensity,ratio,xerr=intensity_error,yerr=ratio_error,fmt='o')

plt.show()

