import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi
analytics = 1
arg_true=0

def e_eq(m1,m2,Ms,inner,outer,length):
    eq = np.zeros(length)
    c=10065.2 #speed of light AU/(yr/2pi)
    G=1
    i=0
    while i < length:   #tailored for 2 planet systems!!
        a1 = inner[i,1]
        e2 = outer[i,2]
        a2 = outer[i,1]
        ep_c = np.sqrt(1.0 - e2*e2)
        n1_sq = (G*(Ms + m1))/(a1*a1*a1)
        gamma = (4*a1*a1*n1_sq/(c*c))*(Ms/m2)*((a2/a1)**3)
        num = (5.0/4.0)*(a1/a2)*e2/(ep_c*ep_c)
        den = np.sqrt(a1/a2)*(m1/m2)/ep_c
        eq[i] = num/abs(1.0 - den + gamma*ep_c**3)
        i += 1
    return eq

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
k2 = np.zeros(N)
Q = np.zeros(N)
for i in range(0,N):
    header = fos.readline()
    header = header.split(",")
    mp[i] = float(header[0])
    rp[i] = float(header[1])
    P[i] = float(header[2])
    k2[i] = float(header[3])
    Q[i] = float(header[4])
header = fos.readline()
header = header.split(",")
tide_delay = float(header[0])

names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','omega (argument of periapsis)','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','2ea(de/dt) - (da/dt)',"-3.2(mu')aen(sinphi)","2ea(de/dt) - (da/dt) - 3.2(mu')aen(sinphi)",'Period Ratio (P$_{i+1}$/P$_{i}$) - j/(j+1)','Resonance Plot','G/G0 - 1']
colors=['b','g','m','r','c','y']
data = np.loadtxt(fos, delimiter="	")

#Figure out array index where eccentricity has settled (i.e. post migration)
ini = 0
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
    if analytics == 1:
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
            edot_e = -(21/2.)*(k2[i]/Q[i])*GM3a3*R5a5/mp[i]
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
    p=data[0::2]
    q=data[1::2]
    arr_len = len(p[:,0])
    eq = e_eq(mp[0],mp[1],Ms,p,q,arr_len)
    plt.plot(p[0:arg3,0], eq[0:arg3], 'o', color = 'orange', markeredgecolor='none', ms = 2)

#plt.xlim([p[arg4,0],p[arg3,0]])
#plt.xlim([0,5000000])
plt.ylim([0,0.1])
#plt.ylim([0,1500])
plt.title(''+name)
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
#plt.legend(loc='upper left',prop={'size':10})
plt.show()
