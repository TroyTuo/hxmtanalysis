from function import *
import sys 

#read two file data
print "################## step one ####################"

filename = "select.fits"
hdulist = pf.open(filename)
tb = hdulist[1].data
tdb = tb.field(1) #!!! THIS IS TDB column of data !!! 
hdulist.close()
raw_data = tdb


MJD = 57946.00529042
data = tdb - tdb[0]

print "################## step two ####################"
# setup import parameters

f0 = 6.581773
f1 = -6.73424e-11
f2 = 1.95e-21
fmax = f0 + 0.5e-3
fmin = f0 - 0.5e-3
fstep = 1e-5

print fmin,fmax
plt.figure()
plt.subplot(2,1,1)
p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,0,fstep,errorbar=False,fig=False,bin_cs=20,bin_profile=20)
fbest = f[chi_square.index(max(chi_square))] 
toa = min(raw_data) + ((1/fbest) * p_x[p.index(max(p))])
print 'toa: ',toa
print min(raw_data)
plt.plot(p_x,p)
min_profile_line = [ p_x[x] for x in xrange(len(p)) if p[x]==min(p) ]
for i in min_profile_line:
    plt.axvline(x=i-0.025,color='r',linewidth=0.5)
    plt.axvline(x=i+0.025,color='r',linewidth=0.5)
plt.subplot(2,1,2)
plt.plot(f,chi_square)
plt.show()
