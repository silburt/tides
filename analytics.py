import sys
import matplotlib.pyplot as plt
import numpy as np
import math

#Analyze the analytic solution of the tidal equations to verify it's correctness
rp = 0.01       #Jupiter=0.1, Earth=0.01, Neptune = 0.035
mp = .0000003   #Jupiter=.001, Earth=3e6, Neptune = 5e5
Ms = 1.         #Units of solar mass
a = 0.2         #AU
G = 1.          #AU^3/Msun*yr*2pi
Qp = 10**-3     #k_2/Q
e = 0.1
dt = 10**-2     #yr/2pi
t = 10**2

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
    da_i[x] = -21.*factor
    de_i[x] = -(21./2.)*factor/(e*a)
    a = a + da_i[x]
    e = e + de_i[x]
    a_i[x] = a
    e_i[x] = e
    t_i[x] = t_i[x-1]+dt

avgda = sum(da_i[:])/float(len(da_i[:]))/dt
avgde = sum(de_i[:])/float(len(de_i[:]))/dt

print "Avg. da/dt=%.20f" % avgda
print "Avg. de/dt=%.20f" % avgde
plt.plot(t_i, a_i)
plt.plot(t_i, e_i)
#plt.plot(t_i, da_i)
#plt.plot(t_i, de_i)
plt.ylim([0,0.3])
plt.xlabel('time')
plt.ylabel('d/dt')
plt.show()