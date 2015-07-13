import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
pi = math.pi
analytics = 1
arg_true=0

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
file_name=str(sys.argv[1])
arg1=int(sys.argv[2])
arg2=int(sys.argv[3])
arg3 = -1
arg4 = 0
if len(sys.argv) == 5:
    arg3 = int(sys.argv[4])
elif len(sys.argv) == 6:
    arg3 = int(sys.argv[4])
    arg4 = int(sys.argv[5])

#Get basic system params from header of file
fos = open(''+file_name, 'r')
header = fos.readline()
header = header.split(",")
name=header[0]
Ms = float(header[1])
Rs = float(header[2])
N = int(header[3])
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
    mig[i] = float(header[5])
header = fos.readline()
header = header.split(",")
tide_delay = float(header[0])
#tide_delay=80000

#Load numerical data
#names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','arg. of peri','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','Resonant Angle (phi = 2*X2 - X1 - w1)','Resonant Angle2 (phi2 = 2*X2 - X1 - w2)','Libration Timescale (order of mag.)','Period Ratio (P$_{i+1}$/P$_{i}$) - j/(j+1)','Resonance Plot','G/G0 - 1']
names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','arg. of peri','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','phi1','phi2','phi3','Period Ratio (P$_{i+1}$/P$_{i}$) - j/(j+1)','Resonance Plot','Delta_Yanqin']
colors=['b','g','m','r','c','y']
data = np.loadtxt(fos, delimiter="	")

#Figure out array index where eccentricity has settled (i.e. post migration)
ini = np.amax(mig)
fini = tide_delay
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
    for i in xrange(1,N):   #outer
        p=data[i::N]
        for j in xrange(0,i):   #inner
            q=data[j::N]
            if(abs(p[-1,3]/(q[-1,3]) - 2) < 0.06):
                plt.plot(p[i_tide:-1,arg1], p[i_tide:-1,3]/q[i_tide:-1,3] - 2., 'o'+colors[inc], label='P$_{'+str(i+1)+',ctlg}$ ='+str(round(P[i],2))+' d, P$_{'+str(j+1)+',ctlg}$='+str(round(P[j],2))+' d, m$_{'+str(i+1)+'}$/m$_{'+str(j+1)+'}$='+str(round(mp[i]/mp[j],3)),markeredgecolor='none', markersize=3)
                inc += 1
elif arg2==12:
    q=data[arg1+1::N]
    length=len(q[:,8])
    x = np.zeros(length)
    y = np.zeros(length)
    for j in xrange(0,length):
        R = q[j,2]
        x[j]=R*math.cos(q[j,8])
        y[j]=R*math.sin(q[j,8])
    gradient = q[:,0]/q[-1,0] #normalize time between 0 and 1
    if arg4 < i_tide:
        arg4 = i_tide
    plt.scatter(x[arg4:arg3], y[arg4:arg3], c=gradient[arg4:arg3], cmap=cm.rainbow, lw=0, label='t$_{max}$ = '+str(round(q[-1,0]/1000000.))+' Myr', alpha = 0.7)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.xlabel('e*cos$\phi$')
    plt.ylabel('e*sin$\phi$')
    cbar = plt.colorbar()
    cbar.set_label('t/t$_{max}$')

elif arg2==14:  #full plot of many things
    q=data[arg1+1::N]
    plt.figure(figsize=(14,6))
    gs = gridspec.GridSpec(3,2)
    plt.subplot(gs[0,0])
    plt.scatter(q[i_tide:arg3,0], q[i_tide:arg3,2])
    plt.xlim([0,q[-1,0]])
    plt.title(''+name)
    plt.ylabel('R $\propto$ e2')    #planet 2 eccentricity
    plt.subplot(gs[1,0])
    plt.scatter(q[i_tide:arg3,0], q[i_tide:arg3,8])
    plt.xlim([0,q[-1,0]])
    plt.ylim([-0.2,6.5])
    plt.ylabel('Resonant Angle $\phi$')
    plt.xlabel('time (years)')
    plt.subplot(gs[2,0])
    plt.scatter(q[i_tide:arg3,0], q[i_tide:arg3,3])
    plt.xlim([0,q[-1,0]])
    plt.ylabel('P2 Period (days)')  #planet 2 Period
    plt.xlabel('time (years)')
    plt.subplot(gs[:,1])
    length=len(q[:,8])
    x = np.zeros(length)
    y = np.zeros(length)
    for j in xrange(0,length):
        R = q[j,2]
        x[j]=R*math.cos(q[j,8])
        y[j]=R*math.sin(q[j,8])
    gradient = q[:,0]/q[-1,0] #normalize time between 0 and 1
    if arg4 < i_tide:
        arg4 = i_tide
    plt.scatter(x[arg4:arg3], y[arg4:arg3], c=gradient[arg4:arg3], cmap=cm.rainbow, lw=0, label='t$_{max}$ = '+str(round(q[-1,0]/1000000.))+' Myr', alpha = 0.7)
    plt.axhline(0, color='black')
    plt.axvline(0, color='black')
    plt.xlabel('e*cos$\phi$')
    plt.ylabel('e*sin$\phi$')
    cbar = plt.colorbar()
    cbar.set_label('t/t$_{max}$')

elif arg2==13:  #Delta growth over time
    inc=0
    for i in xrange(1,N):   #outer
        p=data[i::N]
        for j in xrange(0,i):   #inner
            q=data[j::N]
            if(abs(p[-1,3]/(2*q[-1,3]) - 1) < 0.06):
                plt.plot(p[i_tide:-1,arg1], p[i_tide:-1,3]/(2*q[i_tide:-1,3]) - 1., 'o'+colors[inc], label='P$_{'+str(i+1)+',ctlg}$ ='+str(round(P[i],2))+' d, P$_{'+str(j+1)+',ctlg}$='+str(round(P[j],2))+' d, m$_{'+str(i+1)+'}$/m$_{'+str(j+1)+'}$='+str(round(mp[i]/mp[j],3)),markeredgecolor='none', markersize=3)
                inc += 1
                time = np.arange(0,p[-1,0],p[-1,0]/200.)
                a2=np.median(p[i_tide:f_tide,1])
                a1=np.median(q[i_tide:f_tide,1])
                rad = rp[j]*0.00464913 #Solar Radii to AU
                R5a5 = (rad/a1)**5
                GM3a3 = (Ms/a1)**1.5
                edot_e = (9.*pi/2.)*Qp[j]*GM3a3*R5a5/mp[j]
                B=(mp[i]/mp[j])*(a2/a1)**0.5
                D = 36*mp[j]*mp[j]*B*(1+B)*1.4161
                Delta = 4*10**(-5)*time**(1./3.) #Generic
                #Delta = (D*time*edot_e)**(1./3.)    #Lee et al.
                #k2Qterm = Qp[j]**(1./3.)
                #mterm = (mp[j]*1e6*(1./30.))**(1./3.)
                #rterm = (rp[j]*109.17/2.)**(5./3.)
                #Msterm = Ms**(-8./3.)
                #Pterm = (P[j]/5.)**(-13./9.)
                #Bterm = (2*B + 2*B*B)**(1./3.)
                #Delta_mig3 = (0.006*k2Qterm*mterm*rterm*Msterm*Pterm*Bterm)**3 *(time/5e9)
                #Delta_0 = p[i_tide,3]/(2*q[i_tide,3]) - 1
                #Delta = (Delta_mig3 + Delta_0**3)**(1./3.)     #Lithwick & Wu
                plt.plot(time[0:arg3], Delta[0:arg3], 'k-.', linewidth=3, label='theoretical e(t)')
#******************************
else:
    for i in range(0,N): #range(0,N) only goes to N-1
        p=data[i::N]
        plt.plot(p[arg4:arg3,arg1], p[arg4:arg3,arg2], 'o'+colors[i], label='m$_{'+str(i+1)+'}$='+str(round(100*mp[i]/(3*10**(-6)))/100.)+' m$_{earth}$', markeredgecolor='none')
        if arg2==3:
            plt.plot([p[arg4,0],p[arg3,0]], [P[i], P[i]], label='P$_{catalog}$', color='black', linewidth=2)
    mig_fac = 1.25
    max_mig = max(mig)
    if analytics == 1:
        plt.plot([max_mig, max_mig], [min(data[arg4:arg3,arg2]),max(data[arg4:arg3,arg2])], label='migration ends now!', color='red', linewidth=2)
        plt.plot([tide_delay, tide_delay], [min(data[:,arg2]),max(data[:,arg2])], label='tides turned on now!', color='red', linestyle='dashed', linewidth=2)


#Analytics - plot tidal e - assumes that 'a' is constant, which to first order is true.
if arg2==2 and analytics==1:
    time = np.arange(0,p[-1,0],p[-1,0]/200.)
    for i in range(0,N):
        p=data[i::N]
        a=np.median(p[i_tide:f_tide,1])
        rad = rp[i]*0.00464913 #Solar Radii to AU
        R5a5 = (rad/a)**5
        GM3a3 = (Ms/a)**1.5
        if mp[i] > 0.0:
            edot_e = -(9.*pi/2.)*Qp[i]*GM3a3*R5a5/mp[i]
            t_delay = tide_delay     #delay time before tidal exponential decay starts
            e_t = np.e**(edot_e*(time - t_delay))*np.median(p[i_tide:f_tide,2]) #p[i,2] is a constant of integration, initial e.
            if i == 0:
                plt.plot(time[0:arg3], e_t[0:arg3], 'k-.', linewidth=3, label='theoretical e(t)')
            else:
                plt.plot(time[0:arg3], e_t[0:arg3], 'k-.', linewidth=3)
            tau = -1./edot_e
            print 'tau_e(planet '+str(i+1)+') = '+str(round(tau/1000000.,0))+' Myr'
        #print 'tau_e(planet '+str(i+1)+') = '+str(round(tau,0))+' Years'
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
        GM3mp = (1/mp[i])*(Ms*Ms*Ms)**0.5
        Rp5 = rad**5
        term = -(117.*pi/2.)*Qp[i]*GM3mp*Rp5*e*e
        a_t = (term*time + a0)**(2./13.)
        if i==1:
            plt.plot(time, a_t, 'k-.', linewidth=3, label='theoretical a(t)')
        else:
            plt.plot(time, a_t, 'k-.', linewidth=3)

#plt.ylim([9.9, 10.05])
#plt.ylim([4.9, 5.1])

#plt.ylim([3.95,4.01])
#range=0.05
if arg2==2 and analytics == 1:
    plt.ylim([0.0,0.2])
    plt.title(''+name)
if arg2 != 12 and arg2 != 14:
    #plt.ylim([8.03, 8.12])
    plt.xlim([p[arg4,0],p[arg3,0]])
    plt.xlabel('' + names[arg1])
    plt.ylabel('' + names[arg2])
#plt.legend(loc='upper right',prop={'size':10})
plt.show()
