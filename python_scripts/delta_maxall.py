#The purpose of this macro is to search all round11_Mar26migspeedfac files for each system, and find the run which yields the highest eccentricity for that system. I.e. I am calculating the maximum possible amount of migration due to tides alone.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab
from scipy.stats import itemfreq

arg1='2'
arg2='1'
thresh=0.06
path = '../saved_runs/round11_Mar26migspeedfac/migspeedfac'
pathext = ['15','10','6','4','3','2','1.25','1','0.75','0.5','0.25','0.0625']
N_path = len(pathext)
ext = '_migspeedfac'

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names
N_sys = len(systems)
skip = 0
migfac_index = np.zeros(0)    #catalog the position of the max eccentricity for each run
e0_in = np.zeros(0)
e0_out = np.zeros(0)
a0_in = np.zeros(0)
a0_out = np.zeros(0)
delta = np.zeros(0)
D = np.zeros(0)
i=0
j=0
while i < N_sys:
    print systems[i][1]
    inner = systems[i][24]      #inner planet in res
    outer = systems[i+1][24]    #outer planet in res
    N = systems[i][4]           #number of planets in system
    #Get max(e_i) and a_i
    e0run_in = np.zeros(0)      #array for storing initial e0_in for each run
    e0run_out = np.zeros(0)
    a0run_in = np.zeros(0)
    a0run_out = np.zeros(0)
    j_tracker = np.zeros(0)
    for j in xrange(0,N_path):
        fos = open(path+''+pathext[j]+'/orbits_'+systems[i][1]+''+ext+''+pathext[j]+'.txt','r')
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
        MMRcheck = np.zeros(0)
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
                MMRcheck = np.append(MMRcheck, float(temp_out[3])/float(temp_in[3]) - 2.0)
            elif tt > 70000:
                inMMR_med = np.median(MMRcheck)
                inMMR_min = np.min(MMRcheck)
                inMMR_max = np.max(MMRcheck)
                if abs(inMMR_med) < 0.015 and inMMR_min < 0 and inMMR_max > 0: #ensure planets are in resonance
                    e_in = np.median(e0avg_in)
                    e_out = np.median(e0avg_out)
                    a_in = np.median(a0avg_in)
                    a_out = np.median(a0avg_out)
                    e0run_in = np.append(e0run_in, e_in)
                    e0run_out = np.append(e0run_out, e_out)
                    a0run_in = np.append(a0run_in, a_in)
                    a0run_out = np.append(a0run_out, a_out)
                    j_tracker = np.append(j_tracker, j)     #tracks index of migspeedfac
                exit = 1
            inc_in += N
            inc_out += N
    if len(e0run_in > 0):
        emax_index = int(np.argmax(e0run_in))      #index of emax, regardless of which migspeedfac it belongs to.
        migfac_index = np.append(migfac_index,j_tracker[emax_index]) #info of which migspeedfac it belongs to
        print e0run_in, j_tracker
        e_in = e0run_in[emax_index]
        e_out = e0run_out[emax_index]
        a_in = a0run_in[emax_index]
        a_out = a0run_out[emax_index]
        term = ((a_out/a_in)*math.exp(e_in**2 - e_out**2))**1.5 - 2.0
        D = np.append(D, term)
        delta = np.append(delta, systems[i+1][5]/systems[i][5] - 2.0)      #proximity from 2:1 MMR
    i += 2

#count occurrences of where emax for each system is coming from (i.e. which migspeedfac)
#print migfac_index
freq = itemfreq(migfac_index)
N_diff = freq.shape[0]
for i in xrange(0,N_diff):
    print 'migspeedfac'+pathext[int(freq[i][0])]+' yielded max eccentricity '+str(freq[i][1])+' times'

#plot
binwidth = 0.0001
plt.hist(D, color='green', alpha = 0.8, linewidth=1, bins=np.arange(-0.15, 0.25 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{sim,max}$')
plt.hist(delta, color='blue', alpha = 0.8, linewidth=1, bins=np.arange(-0.15, 0.25 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$')
plt.hist(delta - D, color='black', linewidth=2, bins=np.arange(-0.15, 0.25 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$ - $\Delta_{sim,max}$ ($\Delta$-boost req.)')
plt.plot([0,0],[0,1.25], 'r--', linewidth=2, label='Explainable Via Tides')

plt.ylim([0,1.25])
plt.xlabel('$\Delta$', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.title('Max Evolution Due to Tides ($e_{max}$ all ALL migspeedfacs)')
plt.legend(loc='upper left',prop={'size':10})
pylab.savefig(path+'_delta_max.png')
plt.show()

