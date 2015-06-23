#This macro is to compare my exponential growth to other works (Lee et al. 2013; Lithwick & Wu 2012).

import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
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

#planet values
P1 = 5
P2 = 10
Qp1 = 1./1000
mp1 = 1e-4
mp2 = 1e-4
rp1 = 0.03
rp2 = 0.03
e1 = 0.1

#Star values
Ms = 1
Rs = 1

# time
T = 1e10
Npoints = 1e5

#Lee et al. Delta
time = np.arange(0,T,T/Npoints)
a2=calca(P2,Ms,Rs)
a1=calca(P1,Ms,Rs)
rad = rp1*0.00464913 #Solar Radii to AU
R5a5 = (rad/a1)**5
GM3a3 = (Ms/a1)**1.5
edot_e = (21*pi)*Qp1*GM3a3*R5a5/mp1
B=(mp2/mp1)*(a2/a1)**0.5
D = 36*mp1*mp1*B*(1+B)*1.4161
Delta_L = (D*time*edot_e)**(1./3.)    #Lee et al.
plt.plot(time, Delta_L, 'k-.', linewidth=3, label='Lee et al. (2013)')

#Lithwick Wu
f1 = -1.19
adot_a = 2*e1*edot_e
Gamma = (2+B)*(f1*f1*edot_e*B) + adot_a*f1*f1*B*B/2.
Delta_Li = (9/4 * mp1*mp1*Gamma*time + P2/(2*P1) - 1)**(1./3.)
plt.plot(time, Delta_Li, 'k--', linewidth=3, label='Lithwick & Wu (2012)')

print 1/edot_e/1e6

#Silburt et al Delta
e_t = np.e**(edot_e*(-1)*(time))*e1
Delta_S = ((a2/a1)*np.e**(-e_t*e_t))**1.5 - 2
plt.plot(time, Delta_S - Delta_S[0], 'k', linewidth=3, label='This Work')

plt.ylabel('Planet Spacing $\Delta$', fontsize = 13)
plt.xlabel('Time (years)', fontsize = 13)
plt.legend(loc='upper left',prop={'size':10})
plt.show()






