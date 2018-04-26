from fun import ccf
import matplotlib.pyplot as plt
import numpy as np

def peak_search(profile1, profile2, fig_flag=False):
    '''profile1 is the data of observation
    while profile2 is the standard probe profile'''
    if type(profile1) is str:
        profile1 = np.loadtxt(profile1)
        profile1 = (profile1 - min(profile1))/(max(profile1)-min(profile1))
    
    if type(profile2) is str:
        profile2 = np.loadtxt(profile2)
        profile2 = (profile2 - min(profile2))/(max(profile2)-min(profile2))
        y, delay = ccf(profile1, profile2)
    else:
        profile2 = (profile2 - min(profile2))/(max(profile2)-min(profile2))
        y, delay = ccf(profile1, profile2)
    index_max = delay + np.where(profile2 == max(profile2))
    if fig_flag:
        plt.figure()
        plt.plot(profile1,'b')
        plt.plot(np.roll(profile2,delay),'r')

    index_max = np.asarray(index_max)
    return index_max
