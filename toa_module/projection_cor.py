#!/usr/bin/env python
from __future__ import division
from numpy import *
import math
import matplotlib.pyplot as plt
import pyfits as pf
import glob 
def projection_cor(pathlist,fig_flag=False):
    
    ra = 83.63322083 #degrees
    dec = 22.01446111111 #degrees
    ra = math.radians(ra)
    dec = math.radians(dec)
    angL = math.radians(5.7)
    angS = math.radians(1.1)
    
    #direction vector of source
    xj = cos(ra)*cos(dec)
    yj = sin(ra)*cos(dec)
    zj = sin(dec)
    rj = mat([xj, yj, zj]).T
    
    # read attitude file 
    time_qt_all = array([])
    psf3_all = array([])
    for path in pathlist:
        path = path[0:-1]
        if path == '':continue
        print path
        hdulist = pf.open(glob.glob(path + '/../ACS/*Att*')[0])
        tb = hdulist[3].data
        time_qt = tb.field('Time')
        time_qt_all = append(time_qt_all, time_qt)
        Q1 = tb.field('Q1')
        Q2 = tb.field('Q2')
        Q3 = tb.field('Q3')
        hdulist.close()
        Q0 = abs(sqrt(1 - (Q1**2 + Q2**2 + Q3**2)))
    
        #define Rj
        Rj11 = Q0**2 +Q1**2 - Q2**2 - Q3**2
        Rj21 = 2*(Q1*Q2 - Q0*Q3)
        Rj31 = 2*(Q1*Q3 + Q0*Q2)
        Rj12 = 2*(Q1*Q2 + Q0*Q3)
        Rj22 = Q0**2 + Q2**2 -Q1**2 - Q3**2
        Rj32 = 2*(Q2*Q3 - Q0*Q1)
        Rj13 = 2*(Q1*Q3 - Q0*Q2)
        Rj23 = 2*(Q2*Q3 + Q0*Q1)
        Rj33 = Q0**2 + Q3**2 - Q1**2 - Q2**2
    
        # initial parameters
        Rj = []; rb3 = []; 
        alpha = array([]); 
        beta = array([]); 
        theta = array([]); 
        phi = array([]);
        psf3 = array([]);
    
        for i in xrange(len(Q0)):
            Rj.append(mat([[Rj11[i], Rj12[i], Rj13[i]],[Rj21[i], Rj22[i], Rj23[i]],[Rj31[i], Rj32[i], Rj33[i]]]))
    
            # convert J2000 to satellite coordinate
            rb = Rj[i] * rj
            # convertion to three sets of detectors for HE 
            Rb1 = mat([[0, 1, 0],[0, 0, 1],[1, 0, 0]]) #rotating matrix for detectors aligns 0 degree with Z axis
            Rb2 = mat([[0, 0.5, sqrt(3)/2],[0, -sqrt(3)/2, 0.5],[1, 0, 0]]) #rotating matrix for detectors aligns 60 degrees with Z axis
            Rb3 = mat([[0, 0.5, -sqrt(3)/2],[0, sqrt(3)/2, 0.5],[1, 0, 0]]) #rotating matrix for detectors aligns -60 degrees with Z axis
        #    rb1 = Rb1 * rb
        #    rb2 = Rb2 * rb
            rb3.append(Rb3 * rb)
    
        # calculate projection angles on X, Y axis of collimator
            x = rb3[i][0]; y = rb3[i][1]; z = rb3[i][2]
            alpha = append(alpha, math.atan((x/z)))
            beta = append(beta, math.atan((y/z)))
            theta = append(theta, math.atan((sqrt(x**2+y**2)/z)))
            phi = append(phi, math.atan((y/x)))
    
        # calculate psf for different detectors
        psf3 = (1 - abs(tan(alpha))/tan(angL)) * (1 - abs(tan(beta))/tan(angS)) / sqrt((tan(alpha))**2 + (tan(beta))**2 +1 )
        psf3_all = append(psf3_all, psf3)
    #    for i in xrange(len(psf3)):
    #        if psf3[i] <= 0:
    #            #print psf3[i],math.degrees(alpha[i]),math.degrees(beta[i])
    #            #print (1 - abs(tan(alpha[i]))/tan(angL))
    #            #print (1 - abs(tan(beta[i]))/tan(angS))
    if fig_flag:
        plt.plot(time_qt_all/86400+55927,psf3_all)
        plt.show()
    return time_qt_all, psf3_all

def projection_cor_le(path,fig_flag=False):
    
    ra = 83.63322083 #degrees
    dec = 22.01446111111 #degrees
    ra = math.radians(ra)
    dec = math.radians(dec)
    angL = math.radians(6)
    angS = math.radians(1.6)
    
    #direction vector of source
    xj = cos(ra)*cos(dec)
    yj = sin(ra)*cos(dec)
    zj = sin(dec)
    rj = mat([xj, yj, zj]).T
    
    # read attitude file 
#pathlist = open(pathlistname,'r')
    time_qt_all_le = array([])
    psf1_all_le = array([])
    psf2_all_le = array([])
    psf3_all_le = array([])
    print path
    hdulist = pf.open(glob.glob(path + '/../ACS/*Att*')[0])
    tb = hdulist[3].data
    time_qt = tb.field('Time')
    time_qt_all_le = append(time_qt_all_le, time_qt)
    Q1 = tb.field('Q1')
    Q2 = tb.field('Q2')
    Q3 = tb.field('Q3')
    hdulist.close()
    Q0 = abs(sqrt(1 - (Q1**2 + Q2**2 + Q3**2)))
    
    #define Rj
    Rj11 = Q0**2 +Q1**2 - Q2**2 - Q3**2
    Rj21 = 2*(Q1*Q2 - Q0*Q3)
    Rj31 = 2*(Q1*Q3 + Q0*Q2)
    Rj12 = 2*(Q1*Q2 + Q0*Q3)
    Rj22 = Q0**2 + Q2**2 -Q1**2 - Q3**2
    Rj32 = 2*(Q2*Q3 - Q0*Q1)
    Rj13 = 2*(Q1*Q3 - Q0*Q2)
    Rj23 = 2*(Q2*Q3 + Q0*Q1)
    Rj33 = Q0**2 + Q3**2 - Q1**2 - Q2**2
    
    # initial parameters
    Rj = []; rb1 = []; rb2 = []; rb3 = []; 
    alpha1= array([]); 
    beta1= array([]); 
    theta1= array([]); 
    phi1= array([]);
    psf1 = array([]);
    alpha2= array([]); 
    beta2= array([]); 
    theta2= array([]); 
    phi2= array([]);
    psf2 = array([]);
    alpha3= array([]); 
    beta3= array([]); 
    theta3= array([]); 
    phi3= array([]);
    psf3 = array([]);
    
    for i in xrange(len(Q0)):
        Rj.append(mat([[Rj11[i], Rj12[i], Rj13[i]],[Rj21[i], Rj22[i], Rj23[i]],[Rj31[i], Rj32[i], Rj33[i]]]))
    
        # convert J2000 to satellite coordinate
        rb = Rj[i] * rj
        # convertion to three sets of detectors for HE 
#        Rb1 = mat([[0, 0, 1],[0, -1, 0],[1, 0, 0]]) #rotating matrix for detectors aligns 0 degree with Z axis
#        Rb2 = mat([[0, -sqrt(3)/2, 0.5],[0, -0.5, -sqrt(3)/2],[1, 0, 0]]) #rotating matrix for detectors aligns 60 degrees with Z axis
        Rb3 = mat([[0, sqrt(3)/2, 0.5],[0, -0.5, sqrt(3)/2],[1, 0, 0]]) #rotating matrix for detectors aligns -60 degrees with Z axis
#        rb1.append(Rb1 * rb)
#        rb2.append(Rb2 * rb)
        rb3.append(Rb3 * rb)
    
#    # calculate projection angles on X, Y axis of collimator
#        x = rb1[i][0]; y = rb1[i][1]; z = rb1[i][2]
#        alpha1 = append(alpha1, math.atan((x/z)))
#        beta1 = append(beta1, math.atan((y/z)))
#        theta1 = append(theta1, math.atan((sqrt(x**2+y**2)/z)))
#        phi1 = append(phi1, math.atan((y/x)))
#    
#    # calculate projection angles on X, Y axis of collimator
#        x = rb2[i][0]; y = rb2[i][1]; z = rb2[i][2]
#        alpha2 = append(alpha2, math.atan((x/z)))
#        beta2 = append(beta2, math.atan((y/z)))
#        theta2 = append(theta2, math.atan((sqrt(x**2+y**2)/z)))
#        phi2 = append(phi2, math.atan((y/x)))
    
    # calculate projection angles on X, Y axis of collimator
        x = rb3[i][0]; y = rb3[i][1]; z = rb3[i][2]
        alpha3 = append(alpha3, math.atan((x/z)))
        beta3 = append(beta3, math.atan((y/z)))
        theta3 = append(theta3, math.atan((sqrt(x**2+y**2)/z)))
        phi3 = append(phi3, math.atan((y/x)))
    
#    # calculate psf for different detectors
#    psf1 = (1 - abs(tan(alpha1))/tan(angL)) * (1 - abs(tan(beta1))/tan(angS)) / sqrt((tan(alpha1))**2 + (tan(beta1))**2 +1 )
#    psf1_all_le = append(psf1_all_le, psf1)
#    # calculate psf for different detectors
#    psf2 = (1 - abs(tan(alpha2))/tan(angL)) * (1 - abs(tan(beta2))/tan(angS)) / sqrt((tan(alpha2))**2 + (tan(beta2))**2 +1 )
#    psf2_all_le = append(psf2_all_le, psf2)
    # calculate psf for different detectors
    psf3 = (1 - abs(tan(alpha3))/tan(angL)) * (1 - abs(tan(beta3))/tan(angS)) / sqrt((tan(alpha3))**2 + (tan(beta3))**2 +1 )
    psf3_all_le = append(psf3_all_le, psf3)

    if fig_flag:
        plt.plot(time_qt_all_le/86400+55927,psf3_all_le)
        plt.show()
#    return time_qt_all_le, psf1_all_le, psf2_all_le, psf3_all_le
    return time_qt_all_le,  psf3_all_le
