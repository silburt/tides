#The purpose of this macro is to calculate the minimum eccentricity required by the inner planet in order to achieve the current planet spacing observed due to tides alone, given that the planets started in 2:1 MMR, and that the planets had T years to migrate to that position


import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import islice
import random
random.seed(10)
pi = math.pi

def calca(P,Ms,Rs):
    G = 6.67*10**-11
    P_SI = P*24.*60.*60.
    if Ms == 0:
        Ms = Rs**1.25
    mass = Ms*1.989e30
    a = (P_SI**2*G*mass/(4.*math.pi**2))**(1./3.)
    a /= 1.496*10**11
    return a

def calcm(mp,rp):
    solar2earthRp = 109.21
    earth2solarMp = 3e-6
    if mp == 0. and rp < 0.04:    #Weiss & Marcy
        mp = 2.69*(rp*solar2earthRp)**(0.93) * earth2solarMp
    if mp == 0. and rp >= 0.04:     #Jupiter density
        r3 = (rp*695800000.)**3
        mp = 1330*r3/2e30
    return mp

def calcQp(Qp,rp):
    if rp > 2*0.009156 and rp < 0.1:
        Qp = 1./(2.2e4)
    elif rp >= 0.1:
        Qp = 1./(5.4e4)
    else:
        Qp = 1./40.
    return Qp

#plot initial things
binwidth = 0.0001
y_lim = 1.05
max = 40
med_e_dyn = 0.3 #from Kepler_ei_col.py
linestyles=['solid','dashed','dotted','dashdot']

plt.plot([1.0,1.0],[0,y_lim], 'r--', linewidth=2)
p = plt.axvspan(1.0, max+5, facecolor='red', alpha=0.4, label = 'Unfeasible Eccentricity')
plt.plot([med_e_dyn,med_e_dyn],[0,y_lim], 'b--', linewidth = 2)
p = plt.axvspan(med_e_dyn, 1.0, facecolor='b', alpha=0.4, label = '$e_{i,50}$ (Stability)')

#parameters
T_arr = [1e9,5e9,10e9]      #length of time a system has to achieve current delta spacing
N_T = len(T_arr)
mig_out = 1.01          #(a_i/a_f)_outer planet - an assumption

arg1='2'
arg2='1'

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names

systems = np.delete(systems,[2,3,14,15,16,17,18,19,34,35,36,37,46,47,62,63,68,69])

N_sys = len(systems)

#main loop
for k in xrange(0,N_T):
    e_i = np.zeros(0)
    i=0
    j=0
    T = T_arr[k]
    while i < N_sys:
        Pin = systems[i][5]         #observed period of inner planet (ini condition for sims)
        Pout = systems[i+1][5]
        delta = Pout/Pin - 2.0      #proximity from 2:1 MMR
        Rs = systems[i][17]
        Ms = systems[i][15]
        if Ms == 0:
            Ms = Rs**1.25
        inner = systems[i][24]      #inner planet in res
        outer = systems[i+1][24]    #outer planet in res
        N = systems[i][4]           #number of planets in system
        #Main section - calc inner planet's initial e.
        alpha = 0.630
        b = 1.190
        c = 0.428
        mp_in = systems[i][19]*3e-6
        mp_out = systems[i+1][19]*3e-6
        rp_in = systems[i][22]
        rp_out = systems[i+1][22]
        mp_in = calcm(mp_in,rp_in)
        mp_out = calcm(mp_out,rp_out)
        Qp_in = 0
        Qp_in = calcQp(Qp_in,rp_in)
        rad_in = rp_in*0.00464913 #Solar Radii to AU
        a_in = calca(Pin,Ms,Rs)
        a_out = calca(Pout,Ms,Rs)
        R5a5 = (rad_in/a_in)**5
        GM3a3 = (Ms/a_in)**1.5
        tau_e = 1./((9.*pi/2.)*Qp_in*GM3a3*R5a5/mp_in)
        X_in = T/tau_e
        if X_in > 20:
            ecc = 0.005 + random.random()/100
            #print systems[i][1], 'X_in = ',X_in,'e = ',ecc
        else:
            num1 = ((delta + 2)/2)**(2./3.)
            ecc = (math.log(mig_out*num1)/(math.exp(2*X_in) - 1))**0.5
            #print systems[i][1], 'e_min = ', ecc
        e_i = np.append(e_i,ecc)
        j += 1
        i += 2
    plt.hist(e_i, color='black', alpha = 0.8, linewidth=2, linestyle=linestyles[k], bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$e_{min}$, T='+str(T/1e9)+' Gyr')

plt.ylim([0,y_lim])
plt.xlim([0.001,max+5])
plt.xscale('log')
plt.xlabel('$e_{i,min}$', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=10)
#plt.title('Minimum Eccentricity Required by Inner Planet Given $\Delta$ and t='+str(T/1e9)+' Gyr')
plt.legend(loc='upper left',prop={'size':11})
plt.show()

