import sys
import matplotlib.pyplot as plt
import numpy as np
import math
pi = math.pi
G = 1.
analytics = 0

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
file_name=str(sys.argv[1])
arg1=int(sys.argv[2])
arg2=int(sys.argv[3])

#Get basic system params
fos = open('orbits_sys_char.txt', 'r')
content = fos.readlines()
name = content[0].rstrip('\n')
fos.seek(0)
next(fos)
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
fos = open('runs/'+file_name, 'r')
data = np.loadtxt(fos, delimiter="	")
if arg2==9:
    for i in range(0,N-1):
        p=data[i::N]
        q=data[i+1::N]
        plt.plot(p[:,arg1], q[:,7]/p[:,7], ''+colors[i], label='P$_{'+str(i+2)+'}$/P$_{'+str(i+1)+'}$, m$_{'+str(i+2)+'}$/m$_{'+str(i+1)+'}$='+str(mp[i+1]/mp[i]))
else:
    for i in range(0,N): #range(0,N) only goes to N-1
        p=data[i::N]
        plt.plot(p[:,arg1], p[:,arg2], ''+colors[i], label='m$_{'+str(i+1)+'}$='+str(mp[i]/(3*10**(-6)))+' m$_{earth}$')

#Analytics - plot e
if arg2==2 and analytics==1:
    time = np.arange(0,p[-1,0],p[-1,0]/200.)
    for i in range(0,N):
        a=data[i,1]
        rad = rp[i]*0.00464913 #Solar Radii to AU
        R5a5 = rad*rad*rad*rad*rad/(a*a*a*a*a)
        GM3a3 = (G*Ms*Ms*Ms/(a*a*a))**0.5
        lne = -(9.*pi/2.)*Qp[i]*GM3a3*R5a5/mp[i]
        e_t = math.e**(lne*time)*data[i,2] #p[i,2] is a constant of integration
        plt.plot(time, e_t, 'k-.', linewidth=3, label='theoretical pl.'+str(i+1))

#Analytics - plot a
if arg2==1 and analytics ==1:
    time = np.arange(0,p[-1,0],p[-1,0]/200.)
    for i in range(0,N):
        e   = data[i,2]
        rad = rp[i]*0.00464913 #Solar Radii to AU
        a0  = data[i,1]**(13./2.)
        GM3mp = (1/mp[i])*(G*Ms*Ms*Ms)**0.5
        Rp5 = rad**5
        term = -(117.*pi/2.)*Qp[i]*GM3mp*Rp5*e*e
        a_t = (term*time + a0)**(2./13.)
        plt.plot(time, a_t, 'k-.', linewidth=3, label='theoretical pl.'+str(i+1))

#de = -dt*(9.*pi*0.5)*Qp*GM3a3*R5a5*e/mp

#plt.ylim([0.,0.11])
#plt.ylim([0.2025,0.2075])
plt.xlim([600000,640000])
plt.title(''+name)
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper right')
plt.show()