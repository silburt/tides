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

#Get basic system params from header of file
fos = open(''+file_name, 'r')
header = fos.readline()
header = header.split(",")
name=header[0]
Ms = float(header[1])
Rs = float(header[2])
N = int(header[3])
tide_delay = float(header[4])
rp = np.zeros(N)
mp = np.zeros(N)
Qp = np.zeros(N)
for i in range(0,N):
    header = fos.readline()
    header = header.split(",")
    mp[i] = float(header[0])
    rp[i] = float(header[1])
    Qp[i] = float(header[2])

#Load numerical data
names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','arg. of peri','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','Resonant Angle (phi = 2*X2 - X1 - w1)','Resonant Angle2 (phi2 = 2*X2 - X1 - w2)','Resonant Angle3 (phi3 = w2 - w1)','Period Ratio (P$_{i+1}$/P$_{i}$)','Delta (frac. dist. to res.)']
colors=['b','g','m','r','c','y']
data = np.loadtxt(fos, delimiter="	")
if arg2==11:
    for i in range(0,N-1):
        p=data[i::N]
        q=data[i+1::N]
        plt.plot(p[:,arg1], q[:,3]/p[:,3], 'o'+colors[i], label='P$_{'+str(i+2)+',ini}$ ='+str(round(p[0,3],2))+' d, P$_{'+str(i+1)+',ini}$='+str(round(q[0,3],2))+' d, m$_{'+str(i+2)+'}$/m$_{'+str(i+1)+'}$='+str(round(mp[i+1]/mp[i],3)))
    if tide_delay > 1.:
        plt.plot([tide_delay, tide_delay], [1.99,2.01], label='tides turned on now!', color='black', linewidth=4)
elif arg2==12:
    for i in range(0,N-1):
        p=data[i::N]
        q=data[i+1::N]
        plt.plot(p[:,arg1], q[:,3]/(2*p[:,3]) - 1, 'o'+colors[i], label='P$_{'+str(i+2)+',ini}$ ='+str(round(p[0,3],2))+' d, P$_{'+str(i+1)+',ini}$='+str(round(q[0,3],2))+' d, m$_{'+str(i+2)+'}$/m$_{'+str(i+1)+'}$='+str(round(mp[i+1]/mp[i],3)), markeredgecolor='none')
    if tide_delay > 1.:
        plt.plot([tide_delay, tide_delay], [1.99,2.01], label='tides turned on now!', color='black', linewidth=4)
else:
    for i in range(0,N): #range(0,N) only goes to N-1
        p=data[i::N]
        plt.plot(p[:,arg1], p[:,arg2], 'o'+colors[i], label='m$_{'+str(i+1)+'}$='+str(mp[i]/(3*10**(-6)))+' m$_{earth}$', markeredgecolor='none')
    if tide_delay > 1.:
        plt.plot([tide_delay, tide_delay], [min(data[:,arg2]),max(data[:,arg2])], label='tides turned on now!', color='black', linewidth=4)

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
#plt.xlim([30000,40000])
plt.title(''+name)
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left')
plt.show()