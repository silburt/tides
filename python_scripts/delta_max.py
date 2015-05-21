#The purpose of this macro is the opposite of emin.py - i.e. if I take the initial conditions of the semi-major axis and eccentricity, I can predict the maximum evolution under the influence of tides. Also option to test Delta equation (compare numerics to theory).

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab

arg1='2'
arg2='1'
thresh=0.06

#path = '../saved_runs/round26_May8Qpfac200_tideF/'
#ext = '_Qpfac200_tideF_migfac0.8'
path = '../saved_runs/round27_May8Qpfac200/'
ext = '_Qpfac200_migfac0.8'

numth_comp = 1        #Compare numerics and theory

output = open('../reso/Kepler_ei.txt','w')

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names
#remove systems - Kepler-11 (unstable when in res), Kepler-223 (double res), 226 (inner planet in between resonance),  31 (double res), 331 (double res), 80 (inner planet between resonance), 85 (unstable),
systems = np.delete(systems,[2,3,14,15,16,17,18,19,34,35,36,37,46,47,62,63,68,69])

N_sys = len(systems)
skip = 0
e0_in = np.zeros(0)
e0_out = np.zeros(0)
a0_in = np.zeros(0)
a0_out = np.zeros(0)
P0_out = np.zeros(0)
delta = np.zeros(0)
D = np.zeros(0)
Dnumarr = np.zeros(0)
dD = np.zeros(0)    #Difference in Delta between numerics & theory, used when numth_comp ==1
dafin = np.zeros(0) #Difference in final 'a' between numerics & theory, used when numth_comp ==1 (inner planet)
dafout = np.zeros(0)#Difference in final 'a' between numerics & theory, used when numth_comp ==1 (outer planet)
N_avg_outer = 1000   #number of lines to average over when getting final eccentricity/semi-major axis
i=0
j=0
while i < N_sys:
    print systems[i][1]
    delta = np.append(delta, systems[i+1][5]/systems[i][5] - 2.0)      #proximity from 2:1 MMR
    inner = systems[i][24]      #inner planet in res
    outer = systems[i+1][24]    #outer planet in res
    N = systems[i][4]           #number of planets in system
    #Get e_i and a_i
    fos = open(path+'orbits_'+systems[i][1]+''+ext+'.txt','r')
    if numth_comp == 1:         #extract final a and e of inner/outer planets to compare.
        lines = fos.readlines()
        e1_in = np.zeros(0)
        e1_out = np.zeros(0)
        a1_in = np.zeros(0)
        a1_out = np.zeros(0)
        P1_out = np.zeros(0)
        for k in xrange(1,N_avg_outer+1):
            tempin = lines[-1-k*N+inner+1]
            temp_in = tempin.split("\t")
            a1_in = np.append(a1_in,float(temp_in[1]))
            e1_in = np.append(e1_in,float(temp_in[2]))
            tempout = lines[-1-k*N+outer+1]
            temp_out = tempout.split("\t")
            a1_out = np.append(a1_out,float(temp_out[1]))
            e1_out = np.append(e1_out,float(temp_out[2]))
            P1_out = np.append(P1_out,float(temp_out[3]))
        a_fin = np.median(a1_in)
        e_fin = np.median(e1_in)
        a_fout = np.median(a1_out)
        e_fout = np.median(e1_out)
        P_fout = np.median(P1_out)
    else:
        lines = []
        for k in xrange(0,50000):
            lines.append(fos.readline())
    inc_in = N+2+inner
    inc_out = N+2+outer
    e0avg_in = np.zeros(0)
    e0avg_out = np.zeros(0)
    a0avg_in = np.zeros(0)
    a0avg_out = np.zeros(0)
    P0avg_out = np.zeros(0)
    exit = 0
    dt = 0
    while exit != 1:    #now get initial values, before tides are turned on.
        tempin = lines[inc_in]
        temp_in = tempin.split("\t")
        tt = float(temp_in[0])
        if tt > 50000 and tt < 80000:
            tempout = lines[inc_out]
            temp_out = tempout.split("\t")
            e0avg_in = np.append(e0avg_in, float(temp_in[2]))
            e0avg_out = np.append(e0avg_out, float(temp_out[2]))
            a0avg_in = np.append(a0avg_in, float(temp_in[1]))
            a0avg_out = np.append(a0avg_out, float(temp_out[1]))
            P0avg_out = np.append(P0avg_out, float(temp_out[3]))
        elif tt > 80000:
            e_in = np.median(e0avg_in)      #initial values, before tides are turned on
            e_out = np.median(e0avg_out)
            a_in = np.median(a0avg_in)
            a_out = np.median(a0avg_out)
            P_out = np.median(P0avg_out)
            output.write(systems[i][1]+','+str(e_in)+','+str(a_in)+','+str(e_out)+','+str(a_out)+','+str(inner)+','+str(outer)+'\n')
            e0_in = np.append(e0_in, e_in)
            e0_out = np.append(e0_out, e_out)
            a0_in = np.append(a0_in, a_in)
            a0_out = np.append(a0_out, a_out)
            P0_out = np.append(P0_out, P_out)
            a_fin_th = a_in*math.exp(-e_in**2 + e_fin**2)   #theoretical final 'a' values (inner planet)
            a_fout_th = a_out*math.exp(-e_out**2 + e_fout**2) #theoretical final 'a' values (outer planet)
            Dth = (a_fout_th/a_fin_th)**1.5 - 2.0           #Period ratio, delta_theory
            D = np.append(D, Dth)
            if numth_comp == 1:
                Dnum = (a_fout/a_fin)**1.5 - 2              #Period ratio, delta_numerical
                Dnumarr = np.append(Dnumarr, Dnum)
                dD = np.append(dD, (Dnum - Dth))
                dafin = np.append(dafin, (a_fin_th - a_fin)/(a_in - a_fin))
                dafout = np.append(dafout, (a_fout_th - a_fout)/(a_out - a_fout))
                #print e_in, e_fin, e_out, e_fout, Dth, Dnum, (Dnum - Dth)/Dnum, (a_fin_th - a_fin)/a_fin, (a_fout_th - a_fout)/a_fout
                #print e_in, e_fin, a_in, a_in - a_fin, a_fin - a_fin_th, '(outer)',e_out, e_fout, a_out, a_out - a_fout, a_fout_th - a_fout
                #print i,(Dnum - Dth), Dnum, Dth, N
                print i,Dnum - delta[i/2], Dnum, delta[i/2]
            exit = 1
        inc_in += N
        inc_out += N
    i += 2

print 'median e_i,in =', np.median(e0_in)

xmin = -0.08
xmax = 0.02
binwidth = 0.0001
#plt.hist(D, color='green', alpha = 0.8, linewidth=1, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{sim,max}$')
#plt.hist(delta, color='blue', alpha = 0.8, linewidth=1, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$')
#plt.hist(delta - D, color='black', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$ - $\Delta_{sim,max}$ ($\Delta$-boost req.)')

plt.hist(dD, color='black', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label='$\Delta_{num} - \Delta_{th}$')
plt.hist(Dnumarr - delta, color='black', linestyle='dashed',linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label='$\Delta_{num} - \Delta_{obs}$')

#plt.hist(dafin, color='yellow', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$(a_{th,in} - a_{num,in})/ a_{num,in}$')
#plt.hist(dafout, color='orange', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$(a_{th,out} - a_{num,out})/ a_{num,out})$')

#plt.plot([0,0],[0,1.25], 'r--', linewidth=2, label='Explainable Via Tides')
plt.plot([0,0],[0,1.25], 'r--', linewidth=2)

plt.xlim([-0.07,0.01])
plt.ylim([0,1.05])
plt.xlabel('$d\Delta$', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.legend(loc='upper left',prop={'size':13})
#plt.title('Max Evolution Due to Tides ('+ext+')')
pylab.savefig(path+'delta_max.png')
plt.show()

