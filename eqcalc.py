#The purpose of this macro is to calculate the "equilibrium eccentricity" (inner planet's eccentricity) of a planet based on Eq. 15 from Batygin (2009). Also, the option is here to calculate level curves as a function of a 3rd varying parameter, e.g. see how the e_inner vs. k2_inner relationship changes if inner planet mass is varied as well.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
#import colorsys
import matplotlib.cm as cm

def safety_first(params):    #check for overstepping any bounds
    if params[2] > params[5]:
        print '!ERROR! a1 > a2, reduce parameter space. Exiting.'
        print 'a1=',params[2], 'a2=',params[5]
        exit()
    if params[0] > 10 or params[0] < 0.1:
        print '!ERROR! M < 0.1 or M > 10. Reduce parameter space. Exiting.'
        print 'M=',params[0]
        exit()
    if params[6] > 0.99:
        print '!ERROR! e > 1, reduce parameter space. Exiting.'
        print 'e=',params[6]
        exit()
    if params[4] > 2:   #in Jupiter radii
        print '!ERROR! Rp_inner > 2R_Jupiter, reduce parameter space. Exiting.'
        print 'Rp=',params[4]
        exit()
    if params[3] > 20 or params[7] > 100:    #in Jupiter masses
        print '!ERROR! mp_inner or mp_outer > 20M_Jupiter, reduce parameter space. Exiting.'
        print 'm1=',params[3], 'm2=',params[7]
        exit()

def masterloop(num_points_param, params, min_param, max_param, vp_index, contours):
    if contours == 1: #default values - D for default
        k2_of_half_e_D, half_e_D, factor_D, span_k2_D = masterloop(1, params, 1, 1, vp_index, 0)
    k2_of_half_e = np.zeros(num_points_param)   #arrays
    half_e = np.zeros(num_points_param)
    span_k2 = np.zeros(num_points_param)
    factor = np.zeros(num_points_param)
    colourwheel = cm.rainbow(np.linspace(1, 0, num_points_param))
    for j in xrange(0,num_points_param):
        factor[j] = j*(max_param-min_param)/num_points_param + min_param   #==1 if vp_index = -1
        params[vp_index] *= factor[j]
        M = params[0]               #star
        a1 = params[2]              #inner planet, convert to meters from AU,
        m1 = params[3]*J2S_M       #inner planet, convert to kg from Jupiter mass
        r1 = params[4]*J2AU_R        #inner planet, convert to meters from Jupiter radius
        a2 = params[5]              #outer planet, convert to meters from AU,
        e2 = params[6]              #outer planet, eccentricity
        m2 = params[7]*J2S_M       #outer planet, convert to kg from Jupiter mass
        safety_first(params)
        params[vp_index] /= factor[j]
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
        for i in xrange(0,num_points_e):
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
        #colors, naming and plotting
        labels = ''
        if float(j)/8 - j/8 == 0:
            labels = param_names[vp_index]+'$_{, new}$ ='+str(params[vp_index]*factor[j])
        #(red,green,blue) = colorsys.hsv_to_rgb(0.85*float(j)/num_points_param, 1, 1)
        if len(k2b) > 0:  #make sure there's actually useable values
            if contours == 0:
                colours = 'black'
                labels = 'default curve'
            else:
                colours = colourwheel[j]
                #colours = (red,green,blue)
            a[0].plot(k2b,eb,color=colours,label=labels)
            maximum = max(eb)
            minimum = min(eb)
            half = 0.5*(maximum - minimum) + minimum
            for k in xrange(0,num_points_e):
                if half <= eb[k]:
                    half_e[j] = half
                    k2_of_half_e[j] = k2b[k]
                    span_k2[j] = (maximum - minimum)*k2_of_half_e[j]
                    break
            a[0].scatter([k2_of_half_e[j]],[half_e[j]],s=10, color='black')
    if contours == 1:
        #this is for the label
        a[0].scatter([k2_of_half_e[j]],[half_e[j]],s=10, color='black',label='k2 of half-max(e$_{in}$)')
        #plot span vs. factor
        gradient = k2_of_half_e
        pc = a[1].scatter(factor, span_k2/span_k2_D[0], c=gradient, cmap=cm.rainbow, lw=0, label='k2 of half e', alpha = 0.9, vmin=min(k2_of_half_e), vmax=max(k2_of_half_e))
        cbar = fig.colorbar(pc)
        cbar.set_label('k2 of half-max(e$_{in}$)')
    return k2_of_half_e, half_e, factor, span_k2

#ini args
name = str(sys.argv[1])
vp_index = int(sys.argv[2])     #1 = M, 2 = R, 3 = a1, etc.

#             M    R     a1     m1   r1    a2   e2   m2  name
data_bank=[(0.85,0.81,0.04106,2.13,1.074,3.39,0.829,16,'WASP-53'),
           (1.07,1.28,0.03908,0.725,1.422,2.441,0.5667,57.3,'WASP-81'),
           (1.22,1.56,0.04275,0.851,1.28,1.189,0.691,15.2,'HAT-P-13'),
           (1.126,1.529,0.080,0.442,1.18,2.39,0.21,3.71,'KELT-6'),
           (1.0,1.000,0.05000,0.050,0.352,1.00,0.700,15.0,'Neptune')]
#                  M        R       a1        m1        r1         a2        e2        m2
param_values = [(0.1,2.0),(1,1),(0.5,1.75),(0.1,10.0),(0.2,5.0),(0.5,2.0),(0.75,1.2),(0.1,20.0)]
param_names = ['M$_*$','R$_*$','a$_{in}$','m$_{in}$','r$_{in}$','a$_{out}$','e$_{out}$','m$_{out}$']
param_units = ['M$_{\odot}$','R$_{\odot}$','AU','M$_J$','R$_J$','AU','','M$_J$']

for i in xrange(0,len(data_bank)):
    if name == data_bank[i][8]:
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
max_e = 0.1

#multiple figures plotting
fig, a = plt.subplots(nrows=1, ncols=2, figsize=(15,8))
fig.subplots_adjust(left=0.075, right=0.94)
varying_name = ''

#*****MAIN LOOP**************
contour = 1
if contour == 1:
    num_points_par = 30
    min_par = param_values[vp_index][0]
    max_par = param_values[vp_index][1]
    varying_name = ', Varying '+str(param_names[vp_index])+'='+str(params[vp_index])+str(param_units[vp_index])
else:
    num_points_par = 30
    min_par = 1
    max_par = 1
k2_of_half_e, half_e, factor, span_k2 = masterloop(num_points_par, params, min_par, max_par, vp_index, contour)

#making figures look prettay
fig.text(0.5, 0.95, name+varying_name, ha='center', va='center', rotation='horizontal', fontsize=20)
maxspank2=max(span_k2)
a[0].set_xlim([0,1.5])
a[0].set_ylim([0,2*max(half_e)])
a[0].set_xlabel('k$_{2,in}$', fontsize=15)
a[0].set_ylabel('e$_{in}$', fontsize=15)
if max_par - min_par != 0:
    a[0].legend(loc='upper right',prop={'size':10})

#plot span and k2_half(e_inner) as a function of the varied parameter.
a[1].set_ylabel('Span_k2$_{norm.}$ = k2$_{half(e)}$(e$_{in,max}$ - e$_{in,min}$) / span_k2$_{default}$', fontsize=15)
a[1].set_xlabel('factor ('+param_names[vp_index]+'$_{,fac}$ )', fontsize=15)
#a[1].set_xscale('log')
#a[1].set_xlim([0,max_par*1.1])
#a[1].set_ylim([0,maxspank2+0.005])

#eq.16 - fit
Batygin_fit = 0
if name == 'HAT-P-13' and Batygin_fit == 1 and vp_index == 1:
    k2_fit = np.zeros(num_points_e)
    e_fit = np.zeros(num_points_e)
    for i in xrange(0,num_points_e):
        k2 = 1.5*i/num_points_e
        e_fit[i] = 0.0334 - 0.0985*k2 + 0.188*k2*k2 - 0.184*k2*k2*k2 + 0.069*k2*k2*k2*k2
        k2_fit[i] = k2
    plt.plot(k2_fit,e_fit, 'r',label='Batygin fit')

plt.show()