#This macro is to plot the "resonance tugging" effect. Goes into the round18_Apr9TESTP5m/ folder and plots the various test cases.

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi
Gyr = 1000000000.
A_kq = 300      #Acceleration factor - this is for plotting purposes, to plot the *true* time

#names = ['orbits_TESTP5m_Qpfac20000_migfac0.5.txt','orbits_TESTP5ma_Qpfac100_migfac0.2.txt','orbits_TESTP5mb_Qpfac300_migfac0.5.txt','orbits_TESTP5mc_Qpfac500_migfac0.5.txt','orbits_TESTP5mm_Qpfac100_ei0.10.txt','orbits_TESTP5mma_Qpfac200_ei0.10.txt','orbits_TESTP5mmb_Qpfac175_ei0.10.txt','orbits_TESTP5mmc_Qpfac250_ei0.10.txt',]
#names = ['orbits_TESTP4i_Qpfac300_migfac0.9.txt', 'orbits_TESTP4i_Qpfac300_ei0.0.txt']
names = ['orbits_TESTP4i_Qpfac300_K300_migfac0.8.txt', 'orbits_TESTP4ii_Qpfac300_K300_ei0.00.txt','orbits_TESTP4isingle_Qpfac300_K300_ei0.12.txt']
Nfiles = len(names)
#colors=['black', 'dimgray', 'darkgray', 'lightgray', 'black', 'dimgray', 'darkgray', 'lightgray']
colors=['black', 'gray', 'black']

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
        if i > 1:   #used to be > 3
            line = '--'
            width = 3
        if j == 1:
            labelstr = 'm$_{in}$='+str(round(100*mp[j-1]/(3*10**(-6)))/100.)+' m$_{earth}$'
        axes[1-j].plot(A_kq*p[:,0]/Gyr, p[:,3], colors[i], markeredgecolor='none', linestyle=line, linewidth=width, label=labelstr)
        if i > 3:
            axes[1-j].plot(A_kq*p[:,0]/Gyr, p[:,3], colors[i], markeredgecolor='none', linestyle='-.', linewidth=width, label=labelstr)
    print 'completed file', i

axes[0].set_xlim([A_kq*p[f_tide,0]/Gyr,A_kq*p[-1,0]/Gyr])
#axes[1].set_ylim([4.94,5.06])
#axes[0].set_ylim([10.0,10.12])
#axes[1].set_ylim([3.95,4.00])
#axes[0].set_ylim([7.92, 8.02])
axes[1].set_ylim([3.88,4.00])
axes[0].set_ylim([7.9, 8.02])

axes[1].set_xlabel('time (Gyr)', fontsize=13)
fig.text(0.05, 0.3, 'Inner Planet $P$ (days)', ha='center', va='center', rotation='vertical', fontsize=13)
fig.text(0.05, 0.7, 'Outer Planet $P$ (days)', ha='center', va='center', rotation='vertical', fontsize=13)
#axes[1].set_ylabel('Inner Planet Period (days)', fontsize=13)
#axes[0].set_ylabel('Outer Planet Period (days)', fontsize=13)
#axes[0].set_ylim([1.99,2.01])

fig.subplots_adjust(hspace=0.075)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#axes[0].legend(loc='upper right',prop={'size':10})
plt.show()
