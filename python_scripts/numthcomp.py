#The purpose of this macro is to compare the theoretical predictions of a_f = a_i*exp(-e_i^2 + e_f^2) to the numerical values. Works on a per-system basis. I.e. one .txt file is input

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab

#path = 'numthcomp/orbits_TESTP10k8_Qpfac7500.txt'
path = 'numthcomp/orbits_Kepler-120_Qpfac1000.txt'
#path = 'numthcomp/orbits_TESTK57rp_Qpfac7500.txt'

inner = 0       #inner planet
outer = 1       #outer planet
N = 2           #number of planets

N4avg = 50      #number of lines to average e and a over.
mig_off = 30000 #when to start getting the baseline e value
tides_on = 50000 #when to end getting the baseline e value

#macro
ef_out = np.zeros(0)
af_out = np.zeros(0)
ef_in = np.zeros(0)
af_in = np.zeros(0)
ei_out = np.zeros(0)
ai_out = np.zeros(0)
ei_in = np.zeros(0)
ai_in = np.zeros(0)

fos = open(path,'r')
lines = fos.readlines()
exit = 0
inc_in = N+1+inner  #get past the header
inc_out = N+1+outer
while exit != 1:
    tempin = lines[inc_in]
    temp_in = tempin.split("\t")
    tt = float(temp_in[0])
    if tt > mig_off and tt < tides_on:
    #if tt > 4000000 and tt < 4050000:
        tempout = lines[inc_out]
        temp_out = tempout.split("\t")
        ei_in = np.append(ei_in, float(temp_in[2]))
        ei_out = np.append(ei_out, float(temp_out[2]))
        ai_in = np.append(ai_in, float(temp_in[1]))
        ai_out = np.append(ai_out, float(temp_out[1]))
    elif tt > tides_on:
    #elif tt > 4050000:
        e_iin = np.mean(ei_in)       #initial values, before tides are turned on
        e_iout = np.mean(ei_out)
        a_iin = np.mean(ai_in)
        a_iout = np.mean(ai_out)
        exit = 1
    inc_in += N
    inc_out += N
for i in xrange(1,N4avg+1):         #final values
    tempin = lines[-1-i*N+inner+1]
    temp_in = tempin.split("\t")
    af_in = np.append(af_in,float(temp_in[1]))
    ef_in = np.append(ef_in,float(temp_in[2]))
    tempout = lines[-1-i*N+outer+1]
    temp_out = tempout.split("\t")
    af_out = np.append(af_out,float(temp_out[1]))
    ef_out = np.append(ef_out,float(temp_out[2]))
a_fout = np.mean(af_out)
e_fout = np.mean(ef_out)
a_fin = np.mean(af_in)
e_fin = np.mean(ef_in)

#calculations
a_fin_th = a_iin*math.exp(-e_iin**2 + e_fin**2)   #theoretical final 'a' values (inner planet)
a_fout_th = a_iout*math.exp(-e_iout**2 + e_fout**2) #theoretical final 'a' values (outer planet)
dafin = (a_fin_th - a_fin)/(a_iin - a_fin)
dafout = (a_fout_th - a_fout)/(a_iout - a_fout)
Dth = (a_fout_th/a_fin_th)**1.5 - 2.0
Dnum = (a_fout/a_fin)**1.5 - 2
dD = (Dnum - Dth)/Dth

print ''
print '***outer planet:****'
print 'a_i (num) = ',a_iout, 'e_i (num) = ',e_iout
print 'af_th =',a_fout_th, 'a_i - af_th =', a_iout - a_fout_th
print 'af_num =',a_fout, 'a_i - af_num =', a_iout - a_fout
print 'error between th. and num. - D(af) =',dafout
print '***inner planet:***'
print 'a_i (num) = ',a_iin, 'e_i (num) = ',e_iin
print 'af_th =',a_fin_th, 'a_i - af_th =', a_iin - a_fin_th
print 'af_num =',a_fin, 'a_i - af_num =', a_iin - a_fin
print 'error between th. and num. - D(af) =',dafin
#print '***Delta (distance from resonance)***'
#print 'D(Theory) =', Dth
#print 'D(Numerical) =', Dnum
#print 'Error in D, dD =',dD

