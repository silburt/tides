#The purpose of this macro is the opposite of emin.py - i.e. if I take the initial conditions of the semi-major axis and eccentricity, I can predict the maximum evolution under the influence of tides.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab

arg1='2'
arg2='1'
thresh=0.06
path = '../saved_runs/round11_Mar26migspeedfac0.0625/'
ext = '_migspeedfac0.0625'

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names
N_sys = len(systems)
skip = 0
e0_in = np.zeros(0)
e0_out = np.zeros(0)
a0_in = np.zeros(0)
a0_out = np.zeros(0)
delta = np.zeros(0)
D = np.zeros(0)
i=0
j=0
while i < N_sys:
    delta = np.append(delta, systems[i+1][5]/systems[i][5] - 2.0)      #proximity from 2:1 MMR
    inner = systems[i][24]      #inner planet in res
    outer = systems[i+1][24]    #outer planet in res
    N = systems[i][4]           #number of planets in system
    #Get e_i and a_i
    fos = open(path+'orbits_'+systems[i][1]+''+ext+'.txt','r')
    lines = []
    for k in range(0, 50000):
        lines.append(fos.readline())
    length = len(lines)
    inc_in = N+1+inner
    inc_out = N+1+outer
    e0avg_in = np.zeros(0)
    e0avg_out = np.zeros(0)
    a0avg_in = np.zeros(0)
    a0avg_out = np.zeros(0)
    exit = 0
    while exit != 1:
        tempin = lines[inc_in]
        temp_in = tempin.split("\t")
        tt = float(temp_in[0])
        if tt > 50000 and tt < 70000:
            tempout = lines[inc_out]
            temp_out = tempout.split("\t")
            e0avg_in = np.append(e0avg_in, float(temp_in[2]))
            e0avg_out = np.append(e0avg_out, float(temp_out[2]))
            a0avg_in = np.append(a0avg_in, float(temp_in[1]))
            a0avg_out = np.append(a0avg_out, float(temp_out[1]))
        elif tt > 70000:
            e_in = np.median(e0avg_in)
            e_out = np.median(e0avg_out)
            a_in = np.median(a0avg_in)
            a_out = np.median(a0avg_out)
            e0_in = np.append(e0_in, e_in)
            e0_out = np.append(e0_out, e_out)
            a0_in = np.append(a0_in, a_in)
            a0_out = np.append(a0_out, a_out)
            term = ((a_out/a_in)*math.exp(e_in**2 - e_out**2))**1.5 - 2.0
            D = np.append(D, term)
            exit = 1
        inc_in += N
        inc_out += N
    i += 2

binwidth = 0.0001
plt.hist(D, color='green', alpha = 0.8, linewidth=1, bins=np.arange(-0.15, 0.25 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{sim,max}$')
plt.hist(delta, color='blue', alpha = 0.8, linewidth=1, bins=np.arange(-0.15, 0.25 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$')
plt.hist(delta - D, color='black', linewidth=2, bins=np.arange(-0.15, 0.25 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$ - $\Delta_{sim,max}$ ($\Delta$-boost req.)')
plt.plot([0,0],[0,1.25], 'r--', linewidth=2, label='Explainable Via Tides')

plt.ylim([0,1.25])
plt.xlabel('$\Delta$', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.title('Max Evolution Due to Tides ('+ext+')')
plt.legend(loc='upper left',prop={'size':10})
pylab.savefig(path+'delta_max.png')
plt.show()

