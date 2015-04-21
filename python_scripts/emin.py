#The purpose of this macro is to calculate the minimum eccentricity required by the inner planet in order to achieve the current planet spacing observed due to tides alone (assuming the planets started in the 2:1 MMR, and assuming that the inner planet moves inward according to the current 20Gyr sims in: saved_runs/round9_Mar16Qpfac10000/). This value is then compared to simulated eccentricity values just after migration is finished (pre-tides).

#The calculation is as follows: currently, my simulations are set up such that the inner planet is placed at the current Kepler observed period value (this might have to change), and the outer planet is migrated inwards into resonance.
    #1) Thus, the initial period of the inner planet = P_obs
#Starting the planets in 2:1 MMR, we want to know what the eccentricity the inner planet have to be in order to migrate inwards (due to tides alone) and create the current observed delta = P_out/P_in - 2. But, since the outer planet is also migrating inwards, the inner planet has to migrate even more to obtain the observed delta value.
    #2) Thus, the final period value of the inner planet = P_out_final/(delta + 2.0), where P_out_final is the final period obtained from 20Gyr simulations (saved_runs/round9_Mar16Qpfac10000/), and delta is known from the current observed P2_obs/P1_obs - 2 ratio.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import islice

def get_line(f, n):
    for line in islice(f, n-1, n):
        data = line
    return data

def calca(P,Ms,Rs):
    G = 6.67*10**-11
    P_SI = P*24.*60.*60.
    if Ms == 0:
        Ms = Rs**1.25
    mass = Ms*1.989e30
    a = (P_SI**2*G*mass/(4.*math.pi**2))**(1./3.)
    a /= 1.496*10**11
    return a

arg1='2'
arg2='1'
thresh=0.06
path = '../saved_runs/round8_Mar16Qpfac1/'
ext = '_Qpfac1'
mean_dPout = 0.5 #0.5*(typical amount (in days) that the outer planet migrates in by)

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names
N_sys = len(systems)
skip = 0
e_i = np.zeros(0)
e0 = np.zeros(0)
i=0
j=0
while i < N_sys:
    Pin = systems[i][5]         #observed period of inner planet (ini condition for sims)
    Pout = systems[i+1][5]
    delta = Pout/Pin - 2.0      #proximity from 2:1 MMR
    Ms = systems[i][15]
    Rs = systems[i][17]
    inner = systems[i][24]      #inner planet in res
    outer = systems[i+1][24]    #outer planet in res
    N = systems[i][4]           #number of planets in system
    #outer planet migrates in, find out how much (over 20 Gyr)
    file = '../saved_runs/round9_Mar16Qpfac10000/orbits_'+systems[i][1]+'_Qpfac10000.txt'
    length = num_lines = sum(1 for line in open(file))
    f = open(file,'r')
    temp = get_line(f, length - N + outer + 1)
    data = temp.split("\t")
    Pout_f = float(data[3])     #observed 'P' of outer planet after 20Gyr
    #Method 1 -
    #ain_i = calca(Pin,Ms,Rs)    #ini position for sim, inner planet (outer mig. in).
    #ain_f = calca(Pout_f/(delta + 2.0),Ms,Rs) #final required position of inner, based on current position of outer after 20Gyr of sim
    #****Method 2**** - Assume outer planet doesn't move at all, inner planet starts at Pout/2.0 and does all the migrating in via tides.
    ain_i = calca((Pout + mean_dPout)/2.0,Ms,Rs)
    ain_f = calca(Pin,Ms,Rs)
    ratio = ain_i/ain_f
    if ratio > 1.0:
        term = (math.log(ratio))**0.5
        e_i = np.append(e_i, term)    #from a'/a = 2ee'
        print systems[i][1],', e_min (inner planet) = ', term
    else:
        print '**Warning:'+systems[i][1]+' has ain_i/ain_f < 1.0, removed from analysis**'
        skip = 1
    #compare to simulated eccentricity when run starts
    fos = open(path+'orbits_'+systems[i][1]+''+ext+'.txt','r')
    lines = []
    for k in range(0, 6000):
        lines.append(fos.readline())
    length = len(lines)
    inc = N+1+inner
    e0avg = np.zeros(0)
    exit = 0
    while exit != 1:
        temp2 = lines[inc]
        tempt = temp2.split("\t")
        tt = float(tempt[0])
        if tt > 70000 and tt < 80000:
            temp2out = lines[inc + (outer - inner)]
            temptout = temp2out.split("\t")
            e0avg = np.append(e0avg, float(tempt[2]))
        elif tt > 80000:
            if skip != 1:
                e0 = np.append(e0, np.median(e0avg))
            exit = 1
        inc += N
    j += 1
    i += 2
    skip = 0

binwidth = 0.0001
max = max(e_i) + 0.02
#plt.hist(e0, color='green', alpha = 0.8, linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{sim}$')
plt.hist(e_i, color='blue', alpha = 0.8, linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{min}$')
#plt.hist(e_i - e0, color='black', linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{min}$ - e$_{sim}$ (e-boost req.)')

plt.ylim([0,1.25])
plt.xlabel('e', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.title('Difference in Required e to Explain Obs Period Due to Tides Alone')
plt.legend(loc='upper left',prop={'size':10})
plt.show()

