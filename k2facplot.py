#This macro is made to plot how changing k2fac for a specific system affects the equilibrium eccentricity evolution over time.

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
import glob
pi = math.pi

def label(filename):
    split = filename.split("_")
    k2fac = split[5]
    return [k2fac[5:],split[3],split[4]]


files = glob.glob("saved_runs/amaury/round1_Aug6k2fac/*TESTP2*.txt")
Nfiles = len(files)  #number of files we're dealing with

colors=['purple','blue','skyblue','green','gold','sandybrown','red','brown','black','hotpink']

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
massr = np.zeros(0)

#Get basic system params from header of file
for i in xrange(0,Nfiles):
    fos = open(files[i], 'r')
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
    for j in range(0,N): #range(0,N) only goes to N-1
        p=data[j::N]
        labelstr=''
        line='-'
        width=1
        if j == 1:
            names = label(files[i])
            labelstr = 'k2fac = '+names[0]
        axes[1-j].plot(p[:,0], p[:,2], colors[i], markeredgecolor='none', linestyle=line, linewidth=width, label=labelstr)
    print 'completed file', i

axes[0].set_xlim([0,p[-1,0]])
axes[1].set_ylim([0,0.1])
axes[0].set_ylim([0.2,0.42])
axes[1].set_xlabel('time (years)', fontsize=13)
fig.text(0.05, 0.3, 'Inner Planet $e$', ha='center', va='center', rotation='vertical', fontsize=13)
fig.text(0.05, 0.7, 'Outer Planet $e$', ha='center', va='center', rotation='vertical', fontsize=13)
axes[0].set_title('System='+names[1]+', '+names[2])
#axes[1].set_ylabel('Inner Planet Period (days)', fontsize=13)
#axes[0].set_ylabel('Outer Planet Period (days)', fontsize=13)
#axes[0].set_ylim([1.99,2.01])

fig.subplots_adjust(hspace=0.075)
plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
axes[0].legend(loc='upper right',prop={'size':10})
plt.show()
