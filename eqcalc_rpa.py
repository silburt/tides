#The purpose of this macro is to calculate how the span and k2_at_half_emax change as a function of (Rp/a)_inner, varying a third parameter sometimes in the process.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
#import colorsys
import matplotlib.cm as cm

def safety_first(params):    #check for overstepping any bounds
    if params[3] > params[6]:
        print '!ERROR! a1 > a2, reduce parameter space. Exiting.'
        print 'a1=',params[3], 'a2=',params[6]
        exit()
    if params[1] > 10 or params[1] < 0.1:
        print '!ERROR! M < 0.1 or M > 10. Reduce parameter space. Exiting.'
        print 'M=',params[1]
        exit()
    if params[7] > 0.99:
        print '!ERROR! e > 1, reduce parameter space. Exiting.'
        print 'e=',params[7]
        exit()
    if params[5] > 2:   #in Jupiter radii
        print '!ERROR! Rp_inner > 2R_Jupiter, reduce parameter space. Exiting.'
        print 'Rp=',params[5]
        exit()
    if params[4] > 20 or params[8] > 20:    #in Jupiter masses
        print '!ERROR! mp_inner or mp_outer > 20M_Jupiter, reduce parameter space. Exiting.'
        print 'm1=',params[4], 'm2=',params[8]
        exit()

def masterloop(num_points_param, params, min_param, max_param, vp_index):
    param_names = ['dummy','M$_*$','R$_*$','a$_{in}$','m$_{in}$','r$_{in}$','a$_{out}$','e$_{out}$','m$_{out}$']
    min_rpa_fac = 0.5
    max_rpa_fac = 2
    num_points_rpa = 30
    rpa_fac = np.zeros(num_points_rpa)
    colourwheel = cm.rainbow(np.linspace(1, 0, num_points_param))
    for j in xrange(0,num_points_rpa):  #do this once
        rpa_fac[j] = j*(max_rpa_fac-min_rpa_fac)/num_points_rpa + min_rpa_fac   #==1 if vp_index = -1
    for k in xrange(0,num_points_param):    #parameter varying, e.g. M*
        factor = k*(max_param-min_param)/num_points_param + min_param
        k2_of_half_e = np.zeros(num_points_rpa)   #arrays
        span_k2 = np.zeros(num_points_rpa)
        rp_a = np.zeros(num_points_rpa)
        for j in xrange(0,num_points_rpa):    #varying Rp/a in this loop
            params[vp_index] *= factor
            M = params[1]              #star
            a1 = params[3]*rpa_fac[j]  #inner planet, convert to meters from AU,
            m1 = params[4]*J2S_M       #inner planet, convert to kg from Jupiter mass
            r1 = params[5]*J2AU_R      #inner planet, convert to meters from Jupiter radius
            a2 = params[6]             #outer planet, convert to meters from AU,
            e2 = params[7]             #outer planet, eccentricity
            m2 = params[8]*J2S_M       #outer planet, convert to kg from Jupiter mass
            safety_first(params)
            params[vp_index] /= factor
            rp_a[j] = r1/a1
            ec2_inv = 1.0/(1.0 - e2*e2)
            alpha = a1/a2
            alpha3 = alpha**3
            ab_3 = a1**3
            nb = np.sqrt(G*(M+m1)/ab_3)
            nb_3 = nb**3
            nc = np.sqrt(G*(M+m2)/(a2**3))
            cosdw = 1   #could also be -1
            term1 = (1 + 4*e2*e2)*ec2_inv
            a5R5 = (a1/r1)**5 #tidal terms
            a_c2 = (a1/c)**2
            eb = np.zeros(0)     #reset values
            k2b = np.zeros(0)
            num_points_e = 500
            for i in xrange(0,num_points_e):    #varying e in this loop and solving for k2
                e1 = max_e*i/num_points_e+1e-5
                e1_2 = e1*e1
                eb2_inv = 1.0/(1.0 - e1_2)
                eb2_inv5 = eb2_inv**(-5)
                f2 = eb2_inv5*(1 + 1.5*e1_2 + 0.125*e1_2*e1_2)
                wb_GR = 3*nb_3*eb2_inv*a_c2
                wb_s = 0.75*nb*(m2/M)*alpha3*ec2_inv**(1.5)*(1 - 1.25*alpha*(e2/e1)*ec2_inv*cosdw)
                wc_s = 0.75*nc*(m1/M)*alpha*alpha*ec2_inv*ec2_inv*(1 - 1.25*alpha*(e1/e2)*term1*cosdw)
                numerator = 2*(wc_s - wb_s - wb_GR)*a5R5
                denominator = 15*(M/m1)*f2*nb + nb_3*ab_3*eb2_inv*eb2_inv/(G*m1)
                k2val = numerator/denominator
                if k2val >= 0 and k2val <=1.5:
                    k2b = np.append(k2b,k2val)
                    eb = np.append(eb,e1)
            if len(k2b) > 0:  #make sure there's actually useable values
                maximum = max(eb)
                minimum = min(eb)
                half = 0.5*(maximum - minimum) + minimum
                for h in xrange(0,num_points_e):
                    if half <= eb[h]:
                        k2_of_half_e[j] = k2b[h]
                        span_k2[j] = (maximum - minimum)*k2_of_half_e[j]
                        break
        plt.plot(rp_a,span_k2,color=colourwheel[k],label=param_names[vp_index]+'$_{,new}$ = '+param_names[vp_index]+'*'+str(factor))

#ini args
name = str(sys.argv[1])
vp_index = int(sys.argv[2])     #1 = M, 2 = R, 3 = a1, etc.

#            name      M    R     a1     m1   r1    a2   e2   m2
data_bank=[('HAT-P-13',1.22,1.56,0.04275,0.851,1.28,1.189,0.691,15.2),
           ('Earth',1.0,1.0,0.05,1./300.,0.1,1.0,0.5,1.0),      #Jupiter outer planet
           ('Jupiter',1.0,1.0,0.05,1.0,1.0,1.0,0.5,10),         #10x Jupiter outer planet
           ('Neptune',1.0,1.0,0.05,0.05,0.352,1.0,0.5,2.0),     #2x Jupiter outer planet
           ('WASP-53',0.85,0.81,0.04106,2.13,1.074,3.39,0.829,16),
           ('WASP-81',1.07,1.28,0.03908,0.725,1.422,2.441,0.5667,57.3),
           ('KELT-6',1.126,1.529,0.080,0.442,1.18,2.39,0.21,3.71)]
#                           M       R      a1        m1        r1         a2        e2        m2
param_values = [(-1,-1),(0.5,1.5),(1,1),(0.5,1.75),(0.1,10.0),(0.5,2.0),(0.5,2.0),(0.75,1.2),(0.1,20.0)]
param_names = ['dummy','M$_*$','R$_*$','a$_{in}$','m$_{in}$','r$_{in}$','a$_{out}$','e$_{out}$','m$_{out}$']
param_units = ['dummy','M$_{\odot}$','R$_{\odot}$','AU','M$_J$','R$_J$','AU','','M$_J$']

for i in xrange(0,len(data_bank)):
    if name == data_bank[i][0]:
        params = list(data_bank[i])
        break

#********************
#conversion factors & constants
#G=1 units
J2S_M = 0.0009543       #mass Jupiter -> mass sun
J2S_R = 0.10045         #radius Jupiter -> radius sun
J2AU_R = 0.0004673195   #radius Jupiter->AU
c=10065.2               #speed of light AU/(yr/2pi)
G=1                     #AU^3/Ms/(yr/2pi)^2

#e1 boundaries/# points
num_points_e = 500
max_e = 0.15

#*****MAIN LOOP**************
num_points_par = 5
min_par = param_values[vp_index][0]
max_par = param_values[vp_index][1]
varying_name = ', Varying '+str(param_names[vp_index])+'='+str(params[vp_index])+str(param_units[vp_index])
masterloop(num_points_par, params, min_par, max_par, vp_index)

#title settings
plt.xlabel('(rp/a)$_{in}$',fontsize=15)
plt.ylabel('k2$_{half(e)}$(e$_{in,max}$ - e$_{in,min}$)',fontsize=15)
plt.legend(loc='upper right',prop={'size':11})
plt.title(name)
plt.show()