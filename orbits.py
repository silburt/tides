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
a=0.057216999985679    #AU - THIS CHANGES OVER TIME!!!!
G = 1.          #AU^3/Msun*(yr/2pi)
Qp = 0.1/10.     #k_2/Q
pi = math.pi

com=Ms
R5a5 = rp*rp*rp*rp*rp/(a*a*a*a*a)
GM3a3 = (G*com*com*com/(a*a*a))**0.5
lne = -(9.*pi/2.)*Qp*GM3a3*R5a5/mp
time = np.arange(0,10**7,10**5)
e_t = math.e**(lne*time/(2.*pi)) -0.9
plt.plot(time, e_t, 'k-.', label='analytics')

e = 0.09999999874855819426
dt = 0.00782463923683634696
de = -dt*(9.*pi*0.5)*Qp*GM3a3*R5a5*e/mp

#check - a small da isn't the cause of the discrepancy
#a2 = 0.057310
#R5a5_2 = rp*rp*rp*rp*rp/(a2*a2*a2*a2*a2)
#GM3a3_2 = (G*com*com*com/(a2*a2*a2))**0.5
#lne_2 = -(9.*pi/2.)*Qp*GM3a3_2*R5a5_2/mp
#e_t_2 = math.e**(lne_2*time/(2.*pi)) -0.9
#plt.plot(time, e_t_2, 'g', label='analytics2')

#plt.ylim([0.057216,0.057218])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.show()