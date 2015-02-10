import sys
import matplotlib.pyplot as plt
import numpy as np
import math
pi = math.pi
G = 1.
analytics = 1
arg3=0
arg4=0
arg_true=0

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
P = np.zeros(N)
Qp = np.zeros(N)
mig = np.zeros(N)
for i in range(0,N):
    header = fos.readline()
    header = header.split(",")
    mp[i] = float(header[0])
    rp[i] = float(header[1])
    P[i] = float(header[2])
    Qp[i] = float(header[3])
    mig[i] = float(header[5]) + float(header[5])/4. #damping time too

#Load numerical data
names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','arg. of peri','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','Resonant Angle (phi = 2*X2 - X1 - w1)','Resonant Angle2 (phi2 = 2*X2 - X1 - w2)','Resonant Angle3 (phi3 = w2 - w1)','Period Ratio (P$_{i+1}$/P$_{i}$) - j/(j+1)','Resonance Plot']
colors=['b','g','m','r','c','y']
data = np.loadtxt(fos, delimiter="	")

#Figure out array index where eccentricity has settled (i.e. post migration)
ini = np.amax(mig)
fini = ini + 60000
i_tide = 0
f_tide = 0
var = 0
p=data[0::N]
for i in xrange(0,len(p[:,0])):
    if p[i,0] > ini and var==0:
        i_tide = i
        var = 1
    if p[i,0] > fini:
        f_tide = i
        break

if f_tide == i_tide:
    f_tide += 1

#Choices - main body
if arg2==11:
    inc=0
    for i in xrange(1,N):
        p=data[i::N]
        for j in xrange(0,i):
            q=data[j::N]
            if(abs(p[-1,3]/(2*q[-1,3]) - 1) < 0.05):
                plt.plot(p[1:-1,arg1], p[1:-1,3]/q[1:-1,3] - 2., 'o'+colors[inc], label='P$_{'+str(i+1)+',ctlg}$ ='+str(round(P[i],2))+' d, P$_{'+str(j+1)+',ctlg}$='+str(round(P[j],2))+' d, m$_{'+str(i+1)+'}$/m$_{'+str(j+1)+'}$='+str(round(mp[i]/mp[j],3)),markeredgecolor='none', markersize=3)
                inc += 1
                emax = np.max(q[i_tide:f_tide,2])
                emin = np.min(q[i_tide:f_tide,2])
                resbreak_max = (3.*emax**2)**1.5 *Ms/(31.353*mp[i])
                resbreak_min = (3.*emin**2)**1.5 *Ms/(31.353*mp[i])
                print 'res_break min, max = '+str(round(resbreak_min,2))+ ', '+str(round(resbreak_max,2))+ ' (requires > 1)'
                if resbreak_min > 1 or resbreak_max > 1:
                    a=q[1,1]
                    rad = rp[j]*0.00464913 #Solar Radii to AU
                    R5a5 = rad*rad*rad*rad*rad/(a*a*a*a*a)
                    GM3a3 = (G*Ms*Ms*Ms/(a*a*a))**0.5
                    lne = (9.*pi/2.)*Qp[j]*GM3a3*R5a5/mp[j]
                    print 'planets leave resonance in tau_e = '+str(round(1/lne/1000000.,1))+' Myr?'
                    #plt.plot([1/lne, 1/lne],[-0.01,0.01], 'k--', label=r'$\tau_e$ = '+str(round(1/lne/1000000.,1))+' Myr', linewidth=2)

                #if tide_delay > 1.:
                    #plt.plot([tide_delay, tide_delay], [1.95,2.05], label='tides turned on now!', color='black', linewidth=2)
                    #plt.plot([mig[1], mig[1]], [1.95,2.05], label='Migration stops!', color='black', linewidth=2, ls='--')
elif arg2==12: #use arg0:1=planet-1 & 2, 1=planet-2 & 3, etc. color coded into 4 time snapshots.
    q=data[arg1+1::N]
    length=len(q[:,8])
    x = np.zeros(length)
    y = np.zeros(length)
    for j in xrange(0,length):
        R=15.874*q[j,2]
        x[j]=R*math.cos(q[j,8])
        y[j]=R*math.sin(q[j,8])
    nplots=4
    block=int(length/nplots)
    for i in range(0,nplots):
        plt.plot(x[i*block:(i+1)*block],y[i*block:(i+1)*block],'o'+colors[i], label='Resonance Plot '+str(i)+'*t$_{max}$/4< t <'+str(i+1)+'*t$_{max}$/4', markersize=8-i)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
else:
    for i in range(0,N): #range(0,N) only goes to N-1
        p=data[i::N]
        plt.plot(p[:,arg1], p[:,arg2], 'o'+colors[i], label='m$_{'+str(i+1)+'}$='+str(round(100*mp[i]/(3*10**(-6)))/100.)+' m$_{earth}$', markeredgecolor='none')
        if arg2==3:
            plt.plot([p[0,0],p[-1,0]], [P[i], P[i]], label='P$_{catalog}$', color='black', linewidth=2)
    if tide_delay > 1.:
        plt.plot([tide_delay, tide_delay], [min(data[:,arg2]),max(data[:,arg2])], label='tides turned on now!', color='black', linewidth=2)


#Analytics - plot tidal e - assumes that 'a' is constant, which to first order is true.
if arg2==2 and analytics==1:
    time = np.arange(0,p[-1,0],p[-1,0]/200.)
    for i in range(0,N):
        p=data[i::N]
        a=p[1,1]
        rad = rp[i]*0.00464913 #Solar Radii to AU
        R5a5 = rad*rad*rad*rad*rad/(a*a*a*a*a)
        GM3a3 = (G*Ms*Ms*Ms/(a*a*a))**0.5
        lne = -(9.*pi/2.)*Qp[i]*GM3a3*R5a5/mp[i]
        e_t = math.e**(lne*time)*np.mean(p[i_tide:f_tide,2]) #p[i,2] is a constant of integration, initial e.
        if i == 1:
            plt.plot(time, e_t, 'k-.', linewidth=3, label='theoretical e(t)')
        else:
            plt.plot(time, e_t, 'k-.', linewidth=3)
        tau = -1./lne
        print 'tau_e(planet '+str(i+1)+') = '+str(round(tau/1000000.,0))+' Myr'
    print 't(simulation)   = '+str(p[-1,0]/1000000.)+' Myr'

#Analytics - plot tidal a - this assumes though that eccentricity is constant, which is
#            totally not true. So in the end this plot is false.
if arg2==1 and analytics == 10000:
    time = np.arange(0,p[-1,0],p[-1,0]/200.)
    for i in range(0,N):
        p=data[i::N]
        e   = p[1,2]
        rad = rp[i]*0.00464913 #Solar Radii to AU
        a0  = p[1,1]**(13./2.)
        GM3mp = (1/mp[i])*(G*Ms*Ms*Ms)**0.5
        Rp5 = rad**5
        term = -(117.*pi/2.)*Qp[i]*GM3mp*Rp5*e*e
        a_t = (term*time + a0)**(2./13.)
        if i==1:
            plt.plot(time, a_t, 'k-.', linewidth=3, label='theoretical a(t)')
        else:
            plt.plot(time, a_t, 'k-.', linewidth=3)

#de = -dt*(9.*pi*0.5)*Qp*GM3a3*R5a5*e/mp

#plt.yscale('log')
#plt.xscale('log')
#plt.ylim([0.,0.11])
#plt.ylim([0,0.08])
#plt.xlim([-0.2,0.2])
plt.title(''+name)
if arg2==12:
    plt.xlabel('15.874e*cos$\phi$')
    plt.ylabel('15.874e*sin$\phi$')
else:
    plt.xlabel('' + names[arg1])
    plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.show()