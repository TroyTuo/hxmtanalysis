from function import fsearch
from function import pf
from function import plt
from function import np
import argparse
import os

tstep = 3600
#t0 = 178430000
hdulist = pf.open('1509_screen.fits')
tb = hdulist[1].data
tmax = max(tb.field(0))
tmin = min(tb.field(0))
hdulist.close()

time_range = np.arange(tmin,tmax,tstep)
chi_square_all = [0] * len(time_range)
profile_all = [0] * len(time_range)
time_all = [0] * len(time_range)
fbest_all = [0] * len(time_range)
mjd_all = [0] * len(time_range)
toa_all = [0]* len(time_range)

#select time
plt.figure()
for timei in xrange(len(time_range)):
    select_text = 'fselect 1509_screen.fits select.fits "time >= '+str(time_range[timei])+' && time <= '+str(time_range[timei]+tstep)+'" clobber=yes'
    print select_text
    os.system(select_text)
#read two file data

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",help="name of input file")
    args = parser.parse_args()
    #filename = args.input
    filename = 'select.fits'
    hdulist = pf.open(filename)
    tb = hdulist[1].data
    raw_data = tb.field(0)
    tdb = tb.field(4)
    if len(tdb)==0:continue
#    if (max(tdb)-min(tdb))<=400:continue
    hdulist.close()

    #select baried data
    print "################## step one ####################"
    #data_selected = [tdb[x] for x in index] 
    
    MJD = 57946.00529042
    #data = data_selected - data_selected[0]
    t0 = raw_data[0]
    data = tdb - t0
    
    print "################## step two ####################"
    f0 = 6.581293
    f0 = 6.5823
    f1 = -6.7695374e-11
    f2 = 1.9587e-21
    fmax = f0 + 1e-3
    fmin = f0 - 1e-3
    fstep = 0.5e-4
    
    print fmin,fmax
#    plt.figure()
    plt.subplot(2,1,1)
    p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,f2,fstep,errorbar=False,fig=False,pannel=False,bin_cs=20,bin_profile=10)
    fbest = f[chi_square.index(max(chi_square))] 
    toa = t0 + ((1/fbest) * p_x[p.index(max(p))])
    print 'toa: ',toa
    print min(raw_data)
    plt.plot(p_x,p)
    min_profile_line = [ p_x[x] for x in xrange(len(p)) if p[x]==min(p) ]
    #for i in min_profile_line:
    #    plt.axvline(x=i-0.025,color='r',linewidth=0.5)
    #    plt.axvline(x=i+0.025,color='r',linewidth=0.5)
    plt.subplot(2,1,2)
    plt.plot(f,chi_square)
    
    time_all[timei] = min(raw_data)
    profile_all[timei] = p
    chi_square_all[timei] = chi_square
    fbest_all[timei] = fbest
    mjd_all[timei] = t0/86400 + 55927
    toa_all[timei] = toa
print '###############################################'
print '###############################################'

print 'profile'
print profile_all
print '\n'
print 'MJD'
print mjd_all
print '\n'
print 'chi_square'
print chi_square_all
print '\n'
print 'TOA'
print toa_all
print '###############################################'
print '###############################################'
plt.figure()
profile_all = [ x for x in profile_all if x is not 0 ]
plt.imshow(profile_all,cmap='hot',interpolation='nearest')
plt.show()
