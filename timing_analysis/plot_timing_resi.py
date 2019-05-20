#!/usr/bin/env python

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

data0 = np.loadtxt('./results/resi_pre.dat')
data1 = np.loadtxt('./results/resi_post1.dat')
data2 = np.loadtxt('./results/resi_post2.dat')

toa = [x[0] for x in data0]
resi = [x[1] for x in data0]
error = [x[2] for x in data0]
toa = toa + ([x[0] for x in data1])
resi = resi + ([x[1] for x in data1])
error = error + ([x[2] for x in data1])
toa = toa + ([x[0] for x in data2])
resi = resi + ([x[1] for x in data2])
error = error + ([x[2] for x in data2])

toa1 = [x for x,y,z in zip(toa,resi,error) if y<=0.2]
resi1 = [y for x,y,z in zip(toa,resi,error) if y<=0.2]
error1 = [z for x,y,z in zip(toa,resi,error) if y<=0.2]

toas = np.asarray(toa1)
residual = np.asarray(resi1)
error = np.asarray(error1)



toas = toas - 58100
plt.plot(toas, residual,'.', marker='s',markersize=4,color='black')
for i in xrange(len(toas)):
    plt.plot([toas[i],toas[i]],[residual[i]-error[i], residual[i]+error[i]],'-',color='black',linewidth=0.5)
plt.xlabel('MJD-58100',fontsize=18)
plt.ylabel('Timing Residual(rotational period)',fontsize=18)
plt.show()
