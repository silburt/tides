#The purpose of this macro is to calculate the maximum eccentricity required by the inner planet in order to achieve the current planet spacing observed due to tides alone (assuming the planets started in the 2:1 MMR). This makes the assumption that the outer planet does not migrate appreciably over this time, and that the increase in P2/P1 ratio is due to the inner planet alone. A reasonable assumption.

#you should compare these eccentricities to your current initial eccentricities (just after migration).

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

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

systems=np.genfromtxt('reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float)) #system names
N_sys = len(systems)
e_i = np.zeros(N_sys/2)
i=0
j=0
while i < N_sys:
    Pin = systems[i][5]
    Pout = systems[i+1][5]
    Ms = systems[i][15]
    Rs = systems[i][17]
    ain_f = calca(Pin,Ms,Rs)
    ain_i = calca(Pout/2.0,Ms,Rs)
    e_i[j] = (math.log(ain_i/ain_f))**0.5
    j += 1
    i += 2

binwidth = 0.01
plt.hist(e_i, color='black', linewidth=2, bins=np.arange(0., 0.06 + binwidth, binwidth), histtype='step',cumulative='true', normed='true')
print e_i

plt.ylim([0,1.25])
plt.xlim([0,0.5])
plt.xlabel('eccentricity')
plt.ylabel('counts, total='+str(N_sys/2))
plt.title('Initial Eccentricities')
plt.show()

