import sys
import matplotlib.pyplot as plt
import numpy as np
import math

#Q (Goldreich & Soter, 1966):
#    Mercury <  190
#    Venus <    17
#   Earth =    13 (MacDonald, 1964)
#    Mars <     26
#    Jupiter ~  (1 to 2)e5
#    Saturn ~   (6 to 7)e4
#    Uranus >   7.2e4
#    Neptune >  7.2e4
#
#Q & k_2 (Yoder, 1995 --> But actually from Solar System Dynamics)
#Values in brackets have been estimated (see pg. 166, Ch. 4.10)
#   Mercury=   (100),  (0.1), apparently 0.5 (Padovan et al.)
#    Venus  =   (100),  0.25
#    Earth  =   12,     0.299
#    Mars   =   86      0.14
#
#k_2 of Giant planets (Gavrilov & Zharkov, 1977)
#    Jupiter=   0.379
#    Saturn =   0.341
#    Uranus =   0.104
#    Neptune=   0.127
#

#Analyze the analytic solution of the tidal equations to verify it's correctness
rp = 0.1       #Jupiter=0.1, Earth=0.01, Neptune = 0.035
mp = 0.001   #Jupiter=.001, Earth=3e-6, Neptune = 5e-5
Ms = 1.         #Units of solar mass
a = 0.05         #AU
G = 1.          #AU^3/Msun*yr*2pi
Qp = 0.4/10**5     #k_2/Q
e = 0.1
dt = 10**5     #yr/2pi
t = 5*10**9     #yr/2pi

iterations = int(round(t/dt))
a_i = np.zeros(iterations)
e_i = np.zeros(iterations)
t_i = np.zeros(iterations)
da_i = np.zeros(iterations)
de_i = np.zeros(iterations)

#rp needs to be in units of AU
rp = rp*0.00464913

for x in range(0,iterations):
    factor = ((G*(Ms**3))**(0.5))*(rp**5)*e*e*Qp*dt/(mp*a**(11./2.))
    da_i[x] = -(63./2.)*factor
    de_i[x] = -(63./4.)*factor/(e*a)
    a = a + da_i[x]
    e = e + de_i[x]
    a_i[x] = a
    e_i[x] = e
    t_i[x] = t_i[x-1]+dt

avgda = sum(da_i[:])/float(len(da_i[:]))/dt
avgde = sum(de_i[:])/float(len(de_i[:]))/dt
outputa = -a_i[0]/avgda
outpute = -e_i[0]/avgde

print "Avg. tau_a=%f" % outputa
print "Avg. tau_e=%f" % outpute
plt.plot(t_i, a_i, label='a(t)')
plt.plot(t_i, e_i, label='e(t)')

plt.ylim([0,0.3])
plt.xlabel('time (years)')
plt.ylabel('X(t)')
plt.legend(loc='upper left')
plt.show()