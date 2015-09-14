#The purpose of this macro is to calculate the (instantaneous) equilibrium eccentricity of a planet, which is eq. 36 from Mardling (2007).

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import colorsys

#ini
name = 'HAT-P-13'        #default system
vp_index = -1           #vary_param_index
if len(sys.argv) >= 2:
    name = str(sys.argv[1])
if len(sys.argv) >= 3:
    vp_index = int(sys.argv[2])     #1 = M, 2 = R, 3 = a1, etc.

#            name      M    R     a1     m1   r1    a2   e2   m2
data_bank=[('WASP-53',0.85,0.81,0.04106,2.13,1.074,3.39,0.829,16),
           ('WASP-81',1.07,1.28,0.03908,0.725,1.422,2.441,0.5667,57.3),
           ('HAT-P-13',1.22,1.56,0.04275,0.851,1.28,1.189,0.691,15.2), #0.4383
           ('KELT-6',1.126,1.529,0.080,0.442,1.18,2.39,0.21,3.71)]
#                 M    R     a1     m1   r1    a2   e2   m2
param_values = [(0.5,2.0),(1,1),(0.5,2.0),(0.5,2.0),(0.5,2.0),(0.5,2.0),(0.5,2.0),(0.5,2.0)]
param_names = ['M','R','a1','m1','r1','a2','e2','m2']

for i in xrange(0,len(data_bank)):
    if name == data_bank[i][0]:
        params = list(data_bank[i])
        break

#********************
#conversion factors & constants
#SI
G=6.67*10**(-11)        #G in SI m^3/kg/s^2
J2kg_M = 1.898*10**27   #jupiter mass -> kg
S2kg_M = 1.989*10**30   #solar mass ->kg
J2m_R = 69911000        #jupiter radius -> meters
AU2m_a = 149597871000   #AU -> meters
c=299792458             #c in m/s

#choosing Q_inner
J2S_R = 0.10045         #radius Jupiter -> radius sun
rp = params[5]*J2S_R        #For Q only - convert to solar radius from Jupiter radius.
if rp > 2*0.009156 and rp < 0.1:
    Q1 = 2.2e4
elif rp >= 0.1:
    Q1 = 5.4e4
else:
    Q1 = 40;

#which additional parameter is being varied

if vp_index == -1:  #don't vary a 3rd parameter
    num_points_param = 1
    min_param = 1
    max_param = 1
    vp_index = 0    #just so that an array index error doesn't occur
else:
    num_points_param = 50
    min_param = param_values[vp_index][0]
    max_param = param_values[vp_index][1]
    varying_name = ', Varying '+str(param_names[vp_index]+'='+str(params[vp_index]))

#e1 boundaries/# points
num_points_e = 500
max_e = 0.1

#loop
for j in xrange(0,num_points_param):
    factor = j*max_param/num_points_param + min_param   #==1 if vp_index = -1
    params[vp_index] *= factor
    M = params[1]*S2kg_M        #star
    a1 = params[3]*AU2m_a       #inner planet, convert to meters from AU,
    m1 = params[4]*J2kg_M       #inner planet, convert to kg from Jupiter mass
    r1 = params[5]*J2m_R        #inner planet, convert to meters from Jupiter radius
    a2 = params[6]*AU2m_a       #outer planet, convert to meters from AU,
    e2 = params[7]              #outer planet, eccentricity
    m2 = params[8]*J2kg_M       #outer planet, convert to kg from Jupiter mass
    params[vp_index] /= factor
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
    eb = np.zeros(num_points_e)     #reset values
    k2b = np.zeros(num_points_e)
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
        k2b[i] = numerator/denominator
        eb[i] = e1
    #colors, naming and plotting
    labels = ''
    if float(j)/10 - j/10 == 0:
        labels = param_names[vp_index]+'_new = m1 *'+str(factor)
    (red,green,blue) = colorsys.hsv_to_rgb(0.85*float(j)/num_points_param, 1, 1)
    if min(k2b) < 1.5:  #make sure there's actually useable values
        plt.plot(k2b,eb,color=(red,green,blue),label=labels)

plt.xlim([0,1.5])
plt.xlabel('k$_{2,b}$', fontsize=15)
plt.ylabel('e$_b$', fontsize=15)
plt.title(name+varying_name)
plt.legend(loc='upper left',prop={'size':10})

#eq.16 - fit
Batygin_fit = 0
if name == 'HAT-P-13' and Batygin_fit == 1:
    k2_fit = np.zeros(num_points_e)
    e_fit = np.zeros(num_points_e)
    for i in xrange(0,num_points_e):
        k2 = 1.5*i/num_points_e
        e_fit[i] = 0.0334 - 0.0985*k2 + 0.188*k2*k2 - 0.184*k2*k2*k2 + 0.069*k2*k2*k2*k2
        k2_fit[i] = k2
    plt.plot(k2_fit,e_fit, 'r',label='Batygin fit')
    plt.legend(loc='upper left',prop={'size':10})

plt.show()