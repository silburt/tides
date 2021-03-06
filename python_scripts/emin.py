#The purpose of this macro is to calculate the minimum eccentricity required by the inner planet in order to achieve the current planet spacing observed due to tides alone (assuming the planets started in the 2:1 MMR. The equation for Delta is rearranged in terms of the initial inner and outer planet's eccentricity. We then use Goldreich & Schlichting's Eq. 42 to get everything in terms of the inner planet's eccentricity, and solve for that. We assume that the final eccentricity of both the inner and outer planet ~ 0.


import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from itertools import islice

def get_line(f, n):
    for line in islice(f, n-1, n):
        data = line
    return data

def calca(P,Ms,Rs):
    G = 6.67*10**-11
    P_SI = P*24.*60.*60.
    if Ms == 0:
        Ms = Rs**1.25
    mass = Ms*1.989e30
    a = (P_SI**2*G*mass/(4.*math.pi**2))**(1./3.)
    a /= 1.496*10**11
    return a

arg1='2'
arg2='1'
thresh=0.06
path = '../saved_runs/round8_Mar16Qpfac1/'
ext = '_Qpfac1'

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names
N_sys = len(systems)
skip = 0
e_i = np.zeros(0)
e0 = np.zeros(0)
i=0
j=0
while i < N_sys:
    Pin = systems[i][5]         #observed period of inner planet (ini condition for sims)
    Pout = systems[i+1][5]
    delta = Pout/Pin - 2.0      #proximity from 2:1 MMR
    Rs = systems[i][17]
    Ms = systems[i][15]
    if Ms == 0:
        Ms = Rs**1.25
    inner = systems[i][24]      #inner planet in res
    outer = systems[i+1][24]    #outer planet in res
    N = systems[i][4]           #number of planets in system
#Main section - calc inner planet's initial e.
    alpha = 0.630
    b = 1.190
    c = 0.428
    mu_in = systems[i][19]*3e-6/Ms
    mu_out = systems[i+1][19]*3e-6/Ms
    rp_in = systems[i][22]
    rp_out = systems[i+1][22]
    solar2earthRp = 109.21
    earth2solarMp = 3e-6
    if mu_in == 0. and rp_in < 0.04:    #Weiss & Marcy
        mu_in = 2.69*(rp_in*solar2earthRp)**(0.93) * earth2solarMp/Ms
    if mu_in == 0. and rp_in >= 0.04:     #Jupiter density
        r3 = (rp_in*695800000.)**3
        mu_in = 1330*r3/2e30/Ms
    if mu_out == 0. and rp_out < 0.04:    #Weiss & Marcy
        mu_out = 2.69*(rp_out*solar2earthRp)**(0.93) * earth2solarMp/Ms
    if mu_out == 0. and rp_out >= 0.04:     #Jupiter density
        r3 = (rp_out*695800000.)**3
        mu_out = 1330*r3/2e30/Ms
    num = math.log(alpha * (delta + 2)**(2./3.))
    den = (1 - (c/(2*alpha*b)*(mu_in/mu_out))**2)
    term = (num/den)**0.5
    if term < 1.0 and term > 0.:
        e_i = np.append(e_i, term)    #from a'/a = 2ee'
        print systems[i][1],', e_min (inner planet) = ', term
    else:
        print '**Warning:'+systems[i][1]+' has e_min =',term, ', removed from analysis**'
        skip = 1
    #compare to simulated eccentricity when run starts
    fos = open(path+'orbits_'+systems[i][1]+''+ext+'.txt','r')
    lines = []
    for k in range(0, 6000):
        lines.append(fos.readline())
    length = len(lines)
    inc = N+1+inner
    e0avg = np.zeros(0)
    exit = 0
    while exit != 1:
        temp2 = lines[inc]
        tempt = temp2.split("\t")
        tt = float(tempt[0])
        if tt > 70000 and tt < 80000:
            temp2out = lines[inc + (outer - inner)]
            temptout = temp2out.split("\t")
            e0avg = np.append(e0avg, float(tempt[2]))
        elif tt > 80000:
            if skip != 1:
                e0 = np.append(e0, np.median(e0avg))
            exit = 1
        inc += N
    j += 1
    i += 2
    skip = 0

binwidth = 0.0001
max = max(e_i) + 0.02
#plt.hist(e0, color='green', alpha = 0.8, linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{sim}$')
plt.hist(e_i, color='blue', alpha = 0.8, linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{min}$')
#plt.hist(e_i - e0, color='black', linewidth=2, bins=np.arange(0., max + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = 'e$_{min}$ - e$_{sim}$ (e-boost req.)')

plt.ylim([0,1.25])
plt.xlabel('e', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.title('Difference in Required e to Explain Obs Period Due to Tides Alone')
plt.legend(loc='upper left',prop={'size':10})
plt.show()

