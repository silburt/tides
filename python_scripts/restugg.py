#This macro is to plot the "resonance tugging" effect. Goes into the round18_Apr9TESTP5m/ folder and plots the various test cases.

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

names = ['orbits_TESTP5m_Qpfac10000_migfac0.1.txt','orbits_TESTP5ma_Qpfac100_migfac0.2.txt','orbits_TESTP5mb_Qpfac100_migfac0.2.txt','orbits_TESTP5mc_Qpfac100_migfac0.2.txt','orbits_TESTP5mm_Qpfac10000.txt','orbits_TESTP5mma_Qpfac100.txt','orbits_TESTP5mmb_Qpfac100.txt','orbits_TESTP5mmc_Qpfac100.txt',]
Nfiles = len(names)
colors=['b','g','m','r','b','g','m','r']

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
massr = np.zeros(0)

#Get basic system params from header of file
for i in xrange(0,Nfiles):
    fos = open('restugg/'+names[i], 'r')
    header = fos.readline()
    header = header.split(",")
    name=header[0]
    Ms = float(header[1])
    Rs = float(header[2])
    N = int(header[3])
    rp = np.zeros(N)
    mp = np.zeros(N)
    P = np.zeros(N)
    Qp = np.zeros(N)
    mig = np.zeros(N)
    for j in range(0,N):
        header = fos.readline()
        header = header.split(",")
        mp[j] = float(header[0])
        rp[j] = float(header[1])
        P[j] = float(header[2])
        Qp[j] = float(header[3])
        mig[j] = float(header[5])
    header = fos.readline()
    header = header.split(",")
    tide_delay = float(header[0])
    data = np.loadtxt(fos, delimiter="	")
    #Figure out array index where eccentricity has settled (i.e. post migration)
    ini = np.amax(mig)
    fini = tide_delay
    i_tide = 0
    f_tide = 0
    var = 0
    p=data[0::N]
    for j in xrange(0,len(p[:,0])):
        if p[j,0] > ini and var==0:
            i_tide = j
            var = 1
        if p[j,0] > fini:
            f_tide = j
            break
    if f_tide == i_tide:
        f_tide += 1
    for j in range(0,N): #range(0,N) only goes to N-1
        p=data[j::N]
        labelstr=''
        line='-'
        width=1
        if i > 3:
            line = '--'
            width = 3
        if j == 1:
            labelstr = 'm$_{in}$='+str(round(100*mp[j-1]/(3*10**(-6)))/100.)+' m$_{earth}$, m$_{out}$='+str(round(100*mp[j]/(3*10**(-6)))/100.)+' m$_{earth}$'
        axes[1-j].plot(p[:,0], p[:,3], colors[i], markeredgecolor='none', linestyle=line, linewidth=width, label=labelstr)
        if i > 3:
            axes[1-j].plot(p[:,0], p[:,3], colors[i], markeredgecolor='none', linestyle='-.', linewidth=width, label=labelstr)
    print 'completed file', i

axes[0].set_xlim([p[f_tide,0],p[-1,0]])
axes[1].set_ylim([4.92,5.06])
axes[0].set_ylim([9.98,10.21])
axes[1].set_xlabel('time (years)')
axes[1].set_ylabel('Inner Planet Period (days)')
axes[0].set_ylabel('Outer Planet Period (days)')
#axes[0].set_ylim([1.99,2.01])

fig.subplots_adjust(hspace=0.075)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
axes[0].legend(loc='upper right',prop={'size':9})
plt.show()
