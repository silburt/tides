#The purpose of this macro is to calculate how the span*k2_at_half_emax change as a function of the ratio of two variables of my choice, varying a third parameter in the process.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
#import colorsys
import matplotlib.cm as cm

def check_for_conflict(vp_index,xvar_choice):
    yesprint = 0
    if xvar_choice == 0:
        if vp_index == 1 or vp_index == 4:
            yesprint = 1
    if xvar_choice == 1:
        if vp_index == 2 or vp_index == 6:
            yesprint = 1
    if xvar_choice == 2:
        if vp_index == 0 or vp_index == 6:
            yesprint = 1
    if xvar_choice == 3:
        if vp_index == 0 or vp_index == 2:
            yesprint = 1
    if xvar_choice == 4:
        if vp_index == 1 or vp_index == 3:
            yesprint = 1
    if yesprint == 1:
        print ''
        print '!ERROR! Parameter being varied twice, check and retry. Choose again. '
        print 'Exiting.'
        exit()

def safety_first(M,a1,m1,r1,a2,e2,m2):    #check for overstepping any bounds
    if a1 > a2:
        print '!ERROR! a1 > a2, reduce parameter space. Exiting.'
        print 'a1=',a1, 'a2=',a2
        exit()
    if M > 10 or M < 0.1:
        print '!ERROR! M < 0.1 or M > 10. Reduce parameter space. Exiting.'
        print 'M=',M
        exit()
    if e2 > 0.99:
        print '!ERROR! e > 1, reduce parameter space. Exiting.'
        print 'e=',e2
        exit()
    if r1 > 2:   #in Jupiter radii
        print '!ERROR! Rp_inner > 2R_Jupiter, reduce parameter space. Exiting.'
        print 'Rp=',r1
        exit()
    if m1 > 100 or m2 > 100:    #in Jupiter masses
        print '!ERROR! mp_inner or mp_outer > 20M_Jupiter, reduce parameter space. Exiting.'
        print 'm1=',m1, 'm2=',m2
        exit()

def calc_params(params,vp_index,xvar_fac,xvar_choice,factor,J2S_M,J2AU_R):
    params[vp_index] *= factor
    M = params[0]              #star
    a1 = params[1]             #inner planet, convert to meters from AU,
    m1 = params[2]*J2S_M       #inner planet, convert to kg from Jupiter mass
    r1 = params[3]*J2AU_R      #inner planet, convert to meters from Jupiter radius
    a2 = params[4]             #outer planet, convert to meters from AU,
    e2 = params[5]             #outer planet, eccentricity
    m2 = params[6]*J2S_M       #outer planet, convert to kg from Jupiter mass
    params[vp_index] /= factor
    if xvar_choice == 0:
        a1 *= xvar_fac
        xratio = a2/a1
    elif xvar_choice == 1:
        m2 *= xvar_fac
        xratio = m2/m1
    elif xvar_choice == 2:
        m2 *= xvar_fac
        xratio = m2/M
    elif xvar_choice == 3:
        m1 *= xvar_fac
        xratio = m1/M
    elif xvar_choice == 4:
        a1 *= xvar_fac
        xratio = a1/r1
    safety_first(M,a1,m1,r1,a2,e2,m2)
    return M,a1,m1,r1,a2,e2,m2,xratio

def masterloop(num_points_param, params, min_param, max_param, vp_index,xvar_choice,min_xvar_fac,max_xvar_fac):
    param_names = ['M$_*$','a$_{in}$','m$_{in}$','r$_{in}$','a$_{out}$','e$_{out}$','m$_{out}$']
    param_units = ['M$_{\odot}$','AU','M$_J$','R$_J$','AU','','M$_J$']
    num_points_xvar = 100
    xvar_fac = np.zeros(num_points_xvar)
    max_e = 0.2
    num_points_e = 500
    colourwheel = cm.rainbow(np.linspace(1, 0, num_points_param))
    for j in xrange(0,num_points_xvar):  #do this once
        xvar_fac[j] = j*(max_xvar_fac-min_xvar_fac)/num_points_xvar + min_xvar_fac   #==1 if vp_index = -1
    for k in xrange(0,num_points_param):    #parameter varying, e.g. M*
        factor = k*(max_param-min_param)/num_points_param + min_param
        k2_of_half_e = np.zeros(num_points_xvar)   #arrays
        span_k2 = np.zeros(num_points_xvar)
        xvar = np.zeros(num_points_xvar)          #x ratio value being varied
        for j in xrange(0,num_points_xvar):    #varying Rp/a in this loop
            M,a1,m1,r1,a2,e2,m2,xvar[j] = calc_params(params,vp_index,xvar_fac[j],xvar_choice,factor,J2S_M,J2AU_R)
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
        plt.plot(xvar,span_k2,color=colourwheel[k],label=param_names[vp_index]+'$_{,new}$ = '+str(round(100*params[vp_index]*factor)/100)+' '+param_units[vp_index])

#ini args
name = str(sys.argv[1])
vp_index = int(sys.argv[2])     #0 = M, 1 = a1, etc.
xvar_choice = int(sys.argv[3])   #0 = a2/a1, 1=m2/m1, 2=m2/M, 3=m1/M, 4=a1/r1
check_for_conflict(vp_index,xvar_choice)

#             M     a1     m1   r1    a2   e2   m2    name
data_bank=[(1.22,0.04275,0.851,1.28,1.189,0.691,15.2,'HAT-P-13'),
           (1.0,0.05,1./300.,0.1,1.0,0.5,5.0,'Earth'),      #5xJupiter outer planet
           (1.0,0.05,1.0,1.0,1.0,0.5,5.0,'Jupiter'),        #5x Jupiter outer planet
           (1.0,0.05,0.05,0.352,1.0,0.5,2.0,'Neptune'),     #2x Jupiter outer planet
           (0.85,0.04106,2.13,1.074,3.39,0.829,16,'WASP-53'),
           (1.07,0.03908,0.725,1.422,2.441,0.5667,57.3,'WASP-81'),
           (1.126,0.080,0.442,1.18,2.39,0.21,3.71,'KELT-6')]
#                    M        a1        m1        r1         a2        e2        m2
param_values = [(0.5,1.5),(0.5,1.75),(0.1,10.0),(0.5,2.0),(0.5,2.0),(0.75,1.2),(0.1,2.0)]
param_names = ['M$_*$','a$_{in}$','m$_{in}$','r$_{in}$','a$_{out}$','e$_{out}$','m$_{out}$']
xvar_names = ['a$_{out}$/a$_{in}$','m$_{out}$/m$_{in}$','m$_{out}$/M$_*$','m$_{in}$/M$_*$','a$_{in}$/rp$_{in}$']

for i in xrange(0,len(data_bank)):
    if name == data_bank[i][7]:
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

#*****MAIN LOOP**************
num_points_par = 6
min_par = param_values[vp_index][0]
max_par = param_values[vp_index][1]
min_xvar_fac = 0.01
max_xvar_fac = 2
#varying_name = ', Varying '+str(param_names[vp_index])+'='+str(params[vp_index])+str(param_units[vp_index])
masterloop(num_points_par,params,min_par,max_par,vp_index,xvar_choice,min_xvar_fac,max_xvar_fac)

#title settings
plt.xlabel(xvar_names[xvar_choice],fontsize=15)
plt.ylabel('k2*Span = k2$_{half(e)}$(e$_{in,max}$ - e$_{in,min}$)',fontsize=15)
plt.legend(loc='lower right',prop={'size':11})
plt.yscale('log')
plt.xscale('log')
plt.title(name)
plt.show()
