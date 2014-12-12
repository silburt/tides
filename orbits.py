import sys
import matplotlib.pyplot as plt
import numpy as np
import math

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
file_name=str(sys.argv[1])
arg1=int(sys.argv[2])
arg2=int(sys.argv[3])

#Get basic system params
fos = open('orbits_sys_char.txt', 'r')
syschar = np.loadtxt(fos, delimiter=",")
Ms = syschar[0,0]
Rs = syschar[0,1]
N = int(syschar[0,2])
rp = np.zeros(N)
mp = np.zeros(N)
Qp = np.zeros(N)
for i in range(0,N):
    mp[i] = syschar[i+1,0]
    rp[i] = syschar[i+1,1]
    Qp[i] = syschar[i+1,2]

#Load numerical data
names=['time (years)','Semi-Major Axis (AU)','Eccentricity','inclination','Long. Asc. Node','arg. of peri','Mean Longitude','Period (Days)','mean anomaly','period ratio']
colors=['b','g','m','r','c','y']
fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter="	")
for i in range(0,N): #range(0,N) only goes to N-1
    p=data[i::N]
    if arg2 == 9:
        plt.plot(p[:,arg1], p[:,7]/p[:,7], ''+colors[i], label='pl.'+str(i+1))
    else:
        plt.plot(p[:,arg1], p[:,arg2], ''+colors[i], label='pl.'+str(i+1))

#Analytics - plot e
if arg2==2:
    pi = math.pi
    G = 1.
    time = np.arange(0,p[-1,0],p[-1,0]/200.)
    for i in range(0,N):
        a=data[i,1]
        rad = rp[i]*0.00464913 #Solar Radii to AU
        R5a5 = rad*rad*rad*rad*rad/(a*a*a*a*a)
        GM3a3 = (G*Ms*Ms*Ms/(a*a*a))**0.5
        lne = -(9.*pi/2.)*Qp[i]*GM3a3*R5a5/mp[i]
        e_t = math.e**(lne*time)*data[i,2] #p[i,2] is a constant of integration
        plt.plot(time, e_t, 'k-.', linewidth=3, label='theoretical pl.'+str(i+1))

#de = -dt*(9.*pi*0.5)*Qp*GM3a3*R5a5*e/mp

#plt.ylim([0.057216,0.057218])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper right')
plt.show()