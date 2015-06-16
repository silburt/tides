import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi
analytics = 1
arg_true=0

def e_eq(a1,a2,m1,m2,e2):
    ep_c = (1 - e2*e2)**0.5
    num = (5./4.)*(a1/a2)*e2/(ep_c*ep_c)
    den = (a1/a2)**0.5 *(m1/m2)/ep_c
    result = num/abs(1.0 - den)
    print 'e_eq of inner planet is:',result
    return result

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

names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','arg. of peri','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','2ea(de/dt) - (da/dt)',"-3.2(mu')aen(sinphi)","2ea(de/dt) - (da/dt) - 3.2(mu')aen(sinphi)",'Period Ratio (P$_{i+1}$/P$_{i}$) - j/(j+1)','Resonance Plot','G/G0 - 1']
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
    for i in xrange(1,N):
        p=data[i::N]
        for j in xrange(0,i):
            q=data[j::N]
            if(abs(p[-1,3]/(2*q[-1,3]) - 1) < 0.06):
                plt.plot(p[i_tide:-1,arg1], p[i_tide:-1,3]/q[i_tide:-1,3] - 2., 'o'+colors[inc], label='P$_{'+str(i+1)+',ctlg}$ ='+str(round(P[i],2))+' d, P$_{'+str(j+1)+',ctlg}$='+str(round(P[j],2))+' d, m$_{'+str(i+1)+'}$/m$_{'+str(j+1)+'}$='+str(round(mp[i]/mp[j],3)),markeredgecolor='none', markersize=3)
                inc += 1
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
    a = np.zeros(N)
    e = np.zeros(N)
    for i in range(0,N):
        p=data[i::N]
        a[i]=np.median(p[i_tide:f_tide,1])
        rad = rp[i]*0.00464913 #Solar Radii to AU
        R5a5 = (rad/a[i])**5
        GM3a3 = (Ms/a[i])**1.5
        if mp[i] > 0.0:
            edot_e = -(9.*pi/2.)*Qp[i]*GM3a3*R5a5/mp[i]
            t_delay = tide_delay     #delay time before tidal exponential decay starts
            e[i] = np.median(p[i_tide:f_tide,2]) #p[i,2] is a constant of integration, initial e.
            e_t = np.e**(edot_e*(time - t_delay))*e[i]
            if i == 0:
                plt.plot(time[0:arg3], e_t[0:arg3], 'k-.', linewidth=3, label='theoretical e(t)')
            else:
                plt.plot(time[0:arg3], e_t[0:arg3], 'k-.', linewidth=3)
            tau = -1./edot_e
            print 'tau_e(planet '+str(i+1)+') = '+str(round(tau/1000000.,0))+' Myr'
        #print 'tau_e(planet '+str(i+1)+') = '+str(round(tau,0))+' Years'
    print 't(simulation)   = '+str(p[-1,0]/1000000.)+' Myr'
    eq = e_eq(a[0],a[1],mp[0],mp[1],e[1])
    plt.plot([time[0],time[arg3]], [eq,eq], label='e$_{eq,inner}$', linewidth=3)

plt.xlim([p[arg4,0],p[arg3,0]])
plt.title(''+name)
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.show()
