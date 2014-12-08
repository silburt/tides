import sys
import matplotlib.pyplot as plt
import numpy as np
import math

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=int(sys.argv[1])
arg2=int(sys.argv[2])
N=int(sys.argv[3])          #Number of planets
file_name=str(sys.argv[4])
names=['time (years)','Semi-Major Axis (AU)','Eccentricity','inclination','Long. Asc. Node','arg. of peri','Mean Longitude','Period (Days)','mean anomaly','period ratio']
colors=['b','g','m','r','c']
fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter="	")

for i in range(0,N): #range(0,N) only goes to N-1
    p=data[i::N]
    if arg2 == 9:
        plt.plot(p[:,arg1], p[:,7]/p[:,7], ''+colors[i])
    else:
        plt.plot(p[:,arg1], p[:,arg2], ''+colors[i])

#Analytics
rp = 0.0183*0.00464913       #Jupiter=0.1, Earth=0.01, Neptune = 0.035
mp = 0.00003   #Jupiter=.001, Earth=3e-6, Neptune = 5e-5
Ms = 1.         #Units of solar mass
a = 0.057217    #AU
G = 1.          #AU^3/Msun*(yr/2pi)
Qp = 0.1/10.     #k_2/Q
pi = math.pi

com=Ms+mp
R5a5 = rp*rp*rp*rp*rp/(a*a*a*a*a)
GM3a3 = (G*com*com*com/(a*a*a))**0.5
lne = -(9.*pi/2.)*Qp*GM3a3*R5a5/mp
time = np.arange(0,7*10**7,10**5)
e_t = math.e**(lne*time/(2.*pi)) -0.9
plt.plot(time, e_t, 'k-.', label='analytics')

#e = -(9.*pi/2.)*Qp*(rp**5)*((G*Ms**3)**0.5)/(a**6.5)/mp
#plt.plot(time, math.e**(e*time), 'm', label='analytics2')

#plt.ylim([0.,0.3])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.show()