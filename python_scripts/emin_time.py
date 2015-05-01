#The purpose of this macro is to calculate the minimum eccentricity required by the inner planet in order to achieve the current planet spacing observed due to tides alone, given that the planets started in 2:1 MMR, and that the planets had T years to migrate to that position


import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import islice
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

#parameters
T = 5e9    #length of time a system has to achieve current delta spacing
mig_out = 1.01  #(a_i/a_f)_outer planet - an assumption

arg1='2'
arg2='1'

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names

#systems = np.delete(systems,[14,15,16,17,18,19,34,35,36,37,58,59,62,63,68,69])

N_sys = len(systems)
e_i = np.zeros(0)
e_dyn = np.zeros(0)
i=0
j=0
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
    solar2earthRp = 109.21
    earth2solarMp = 3e-6
    if mp_in == 0. and rp_in < 0.04:    #Weiss & Marcy
        mp_in = 2.69*(rp_in*solar2earthRp)**(0.93) * earth2solarMp
    if mp_in == 0. and rp_in >= 0.04:     #Jupiter density
        r3 = (rp_in*695800000.)**3
        mp_in = 1330*r3/2e30
    if mp_out == 0. and rp_out < 0.04:    #Weiss & Marcy
        mp_out = 2.69*(rp_out*solar2earthRp)**(0.93) * earth2solarMp
    if mp_out == 0. and rp_out >= 0.04:     #Jupiter density
        r3 = (rp_out*695800000.)**3
        mp_out = 1330*r3/2e30
    Qp_in = 0
    if rp_in > 2*0.009156 and rp_in < 0.1:
        Qp_in = 1./(2.2e4)
    elif rp_in >= 0.1:
        Qp_in = 1./(5.4e4)
    else:
        Qp_in = 1./40.
    rad_in = rp_in*0.00464913 #Solar Radii to AU
    a_in = calca(Pin,Ms,Rs)
    a_out = calca(Pout,Ms,Rs)
    R5a5 = (rad_in/a_in)**5
    GM3a3 = (Ms/a_in)**1.5
    tau_e = 1./((9.*pi/2.)*Qp_in*GM3a3*R5a5/mp_in)
    X_in = T/tau_e
    if X_in > 20:
        ecc = 0.01
        print systems[i][1], 'X_in = ',X_in,'e ~ 0.01'
    else:
        num1 = ((delta + 2)/2)**(2./3.)
        ecc = (math.log(mig_out*num1)/(math.exp(2*X_in) - 1))**0.5
        print systems[i][1], 'e_min = ', ecc
    e_i = np.append(e_i,ecc)
    R_H = 0.5*(a_in + a_out)*((mp_in + mp_out)/(3*Ms))**(1./3.)
    K = (a_out - a_in)/R_H
    K50 = 0.7*np.log10(T/Pin) + 2.87
    print K50, K
    e_dyn = np.append(e_dyn, (K - K50)*0.01)    #max e from stability arg
    j += 1
    i += 2

binwidth = 0.0001
y_lim = 1.25
max = max(e_i) + 0.02
med_e_dyn = np.median(e_dyn)

plt.plot([1.0,1.0],[0,y_lim], 'r--', linewidth=2)
p = plt.axvspan(1.0, max, facecolor='red', alpha=0.4, label = 'Unfeasible Eccentricity')
plt.plot([med_e_dyn,med_e_dyn],[0,y_lim], 'b--', linewidth = 2)
p = plt.axvspan(med_e_dyn, max, facecolor='b', alpha=0.4, label = 'Max Eccentricity (Stability Arg.)')

plt.hist(e_i, color='black', alpha = 0.8, linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{min}$')

plt.ylim([0,y_lim])
plt.xscale('log')
plt.xlabel('e', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.title('Minimum Eccentricity Required by Inner Planet Given $\Delta$ and t='+str(T/1e9)+' Gyr')
plt.legend(loc='upper left',prop={'size':10})
plt.show()

