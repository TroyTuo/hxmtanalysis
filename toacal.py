from function import *
import sys 
import argparse

#read two file data

parser = argparse.ArgumentParser()
filename = "he_screen.fits"
parser.add_argument("-i","--input",help="name of input file")
#parser.add_argument("-f0","--fbest",help="f0",type=float)
#parser.add_argument("-f1","--f1",help="f1",type=float)
#parser.add_argument("-f2","--f2",help="f2",type=float)
#parser.add_argument("-fstep","--fstep",help="fstep",type=float)
args = parser.parse_args()
filename = args.input
hdulist = pf.open(filename)
tb = hdulist[1].data
tdb = tb.field(1)
hdulist.close()
raw_data = tdb

#select baried data
print "################## step one ####################"
#data_selected = [tdb[x] for x in index] 

MJD = 57946.00529042
#data = data_selected - data_selected[0]
data = tdb - tdb[0]

print "################## step two ####################"
f0 = 6.581293
f0 = 6.581783
f0 = 29.63895251
f1 = -6.7695374e-11
f1 = -3.6873e-10
f2 = 1.9587e-21
fmax = f0 + 0.5e-3
fmin = f0 - 0.5e-3
fstep = 0.5e-4
global_bin_cs = 20
global_bin_profile = 20

print fmin,fmax
plt.figure()
plt.subplot(2,1,1)
p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,0,fstep,errorbar=False,fig=False,bin_cs=global_bin_cs,bin_profile=global_bin_profile)
fbest = f[chi_square.index(max(chi_square))] 
toa = min(raw_data) + ((1/fbest) * p_x[p.index(max(p))])
print 'toa: ',toa
print min(raw_data)
plt.plot(p_x,p)
min_profile_line = [ p_x[x] for x in xrange(len(p)) if p[x]==min(p) ]
for i in min_profile_line:
    plt.axvline(x=i-(0.5/global_bin_profile),color='r',linewidth=0.5);print 'phase_off_0: '+i-(0.5/global_bin_profile)
    plt.axvline(x=i+(0.5/global_bin_profile),color='r',linewidth=0.5);print 'phase_oof_1: '+i+(0.5/global_bin_profile)
plt.subplot(2,1,2)
plt.plot(f,chi_square)
plt.show()
