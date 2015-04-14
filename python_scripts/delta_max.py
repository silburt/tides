#The purpose of this macro is the opposite of emin.py - i.e. if I take the initial conditions of the semi-major axis and eccentricity, I can predict the maximum evolution under the influence of tides. Also option to test Delta equation (compare numerics to theory).

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab

arg1='2'
arg2='1'
thresh=0.06
path = '../saved_runs/round13_Apr1Qpfac100/'
ext = '_Qpfac100'
#path = '../saved_runs/round11_Mar26migspeedfac/migspeedfac1.25/'
#ext = '_migspeedfac1.25'

numth_comp = 1        #Compare numerics and theory

output = open('../reso/Kepler_ei.txt','w')

systems=np.genfromtxt('../reso/full/'+arg1+':'+arg2+'_systems_fulldetail.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,int)) #system names
N_sys = len(systems)
skip = 0
e0_in = np.zeros(0)
e0_out = np.zeros(0)
a0_in = np.zeros(0)
a0_out = np.zeros(0)
delta = np.zeros(0)
D = np.zeros(0)
dD = np.zeros(0)    #Difference in Delta between numerics & theory, used when numth_comp ==1
dafin = np.zeros(0) #Difference in final 'a' between numerics & theory, used when numth_comp ==1 (inner planet)
dafout = np.zeros(0)#Difference in final 'a' between numerics & theory, used when numth_comp ==1 (outer planet)
N_avg_outer = 100   #number of lines to average over when getting final eccentricity/semi-major axis
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
        for k in xrange(1,N_avg_outer+1):
            tempin = lines[-1-k*N+inner+1]
            temp_in = tempin.split("\t")
            a1_in = np.append(a1_in,float(temp_in[1]))
            e1_in = np.append(e1_in,float(temp_in[2]))
            tempout = lines[-1-k*N+outer+1]
            temp_out = tempout.split("\t")
            a1_out = np.append(a1_out,float(temp_out[1]))
            e1_out = np.append(e1_out,float(temp_out[2]))
        a_fin = np.median(a1_in)
        e_fin = np.median(e1_in)
        a_fout = np.median(a1_out)
        e_fout = np.median(e1_out)
    else:
        lines = []
        for k in xrange(0,50000):
            lines.append(fos.readline())
    inc_in = N+1+inner
    inc_out = N+1+outer
    e0avg_in = np.zeros(0)
    e0avg_out = np.zeros(0)
    a0avg_in = np.zeros(0)
    a0avg_out = np.zeros(0)
    exit = 0
    dt = 0
    while exit != 1:    #now get initial values, before tides are turned on.
        tempin = lines[inc_in]
        temp_in = tempin.split("\t")
        tt = float(temp_in[0])
        if tt > 30000 and tt < 55000:
            tempout = lines[inc_out]
            temp_out = tempout.split("\t")
            e0avg_in = np.append(e0avg_in, float(temp_in[2]))
            e0avg_out = np.append(e0avg_out, float(temp_out[2]))
            a0avg_in = np.append(a0avg_in, float(temp_in[1]))
            a0avg_out = np.append(a0avg_out, float(temp_out[1]))
        elif tt > 55000:
            e_in = np.median(e0avg_in)      #initial values, before tides are turned on
            e_out = np.median(e0avg_out)
            a_in = np.median(a0avg_in)
            a_out = np.median(a0avg_out)
            output.write(systems[i][1]+','+str(e_in)+','+str(a_in)+','+str(e_out)+','+str(a_out)+','+str(inner)+','+str(outer)+'\n')
            e0_in = np.append(e0_in, e_in)
            e0_out = np.append(e0_out, e_out)
            a0_in = np.append(a0_in, a_in)
            a0_out = np.append(a0_out, a_out)
            a_fin_th = a_in*math.exp(-e_in**2 + e_fin**2)   #theoretical final 'a' values (inner planet)
            a_fout_th = a_out*math.exp(-e_out**2 + e_fout**2) #theoretical final 'a' values (outer planet)
            Dth = (a_fout_th/a_fin_th)**1.5 - 2.0           #Period ratio, delta_theory
            D = np.append(D, Dth)
            if numth_comp == 1:
                Dnum = (a_fout/a_fin)**1.5 - 2              #Period ratio, delta_numerical
                dD = np.append(dD, (Dnum - Dth)/Dth)
                dafin = np.append(dafin, (a_fin_th - a_fin)/(a_in - a_fin))
                dafout = np.append(dafout, (a_fout_th - a_fout)/(a_out - a_fout))
                #print e_in, e_fin, e_out, e_fout, Dth, Dnum, (Dnum - Dth)/Dnum, (a_fin_th - a_fin)/a_fin, (a_fout_th - a_fout)/a_fout
                #print e_in, e_fin, a_in, a_in - a_fin, a_fin - a_fin_th, '(outer)',e_out, e_fout, a_out, a_out - a_fout, a_fout_th - a_fout
                print (Dnum - Dth)/Dth, Dnum, Dth, N
            exit = 1
        inc_in += N
        inc_out += N
    i += 2

xmin = -5
xmax = 5
binwidth = 0.0001
plt.hist(D, color='green', alpha = 0.8, linewidth=1, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{sim,max}$')
plt.hist(delta, color='blue', alpha = 0.8, linewidth=1, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$')
plt.hist(delta - D, color='black', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$\Delta_{obs}$ - $\Delta_{sim,max}$ ($\Delta$-boost req.)')
plt.hist(dD, color='purple', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$(\Delta_{th} - \Delta_{num})/ \Delta_{num}$ (Does theory agree with num)')
#plt.hist(dafin, color='yellow', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$(a_{th,in} - a_{num,in})/ a_{num,in}$')
#plt.hist(dafout, color='orange', linewidth=2, bins=np.arange(xmin, xmax + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label = '$(a_{th,out} - a_{num,out})/ a_{num,out})$')

plt.plot([0,0],[0,1.25], 'r--', linewidth=2, label='Explainable Via Tides')


plt.ylim([0,1.25])
plt.xlabel('$\Delta$', fontsize=16)
plt.ylabel('cdf, counts='+str(N_sys/2), fontsize=16)
plt.title('Max Evolution Due to Tides ('+ext+')')
plt.legend(loc='upper left',prop={'size':8})
pylab.savefig(path+'delta_max.png')
plt.show()

