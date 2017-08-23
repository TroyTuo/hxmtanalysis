from function import *
import sys 

filename = "/Users/tuoyouli/Desktop/fermi_toa/data/bary_1deg.fits"
hdulist = pf.open(filename)
tb = hdulist[1].data
raw_data = tb.field(9)
hdulist.close() 


print raw_data
MJD = 57946.00529042
data_tmp = raw_data-raw_data[0]

f0 = 29.640491#29.640491 
#f0 = 12.3448
f1=-3.68754e-10
f1 = -3.688111035100001e-10
f2=2.1e-20
fmax = f0 + 1e-5
fmin = f0 - 1e-5


data = data_tmp[0:300]
print fmin,fmax
p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,f2,1e-7,errorbar=True,fig=True)
fbest = f[chi_square.index(max(chi_square))] 
toa = min(raw_data) + ((1/fbest) * p_x[p.index(max(p))])
print 'toa: ',toa
#plt.figure('total')
#plt.subplot(2,1,1)
#plt.plot(p_x,p)
#plt.subplot(2,1,2)
#plt.plot(f,chi_square)
#plt.show()
#

profile_list = []
#for i in range(0,9):
#    tmin = 200*i
#    tmax = 200*(i+1)
#    data = data_tmp[tmin:tmax]
#    print fmin,fmax
#    p_x,p,f,chi_square = fsearch(data,fmin,fmax,f1,f2,1e-7)
#    fbest = f[chi_square.index(max(chi_square))] 
#    toa = min(raw_data) + ((1/fbest) * p_x[p.index(max(p))])
#    print 'toa: ',toa
#    profile_list.append(p)
#    #data1 = raw_data[0:500]
#    #p_x1,p1,f_cut,chi_square_cut = fsearch(data1,fmin-0.0001,fmax+0.0001,f1,f2,1e-8)
#    #fbest = f[chi_square.index(max(chi_square))] 
#    #toa = min(data1) + ((1/fbest) * p_x1[p1.index(max(p1))])
#    #print 'toa: ',toa
#
#    plt.figure('total')
#    plt.subplot(2,1,1)
#    plt.plot(p_x,p)
#    plt.figure(2)
#    plt.subplot(2,1,2)
#    plt.plot(f,chi_square)
#    plt.show()


#plt.figure('cut')
#plt.subplot(2,1,1)
#plt.plot(p_x,p)
#plt.subplot(2,1,2)
#plt.plot(f_cut,chi_square_cut)
   # plt.show()
