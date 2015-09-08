#The purpose of this macro is to calculate the (instantaneous) equilibrium eccentricity of a planet, which is eq. 36 from Mardling (2007).

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

#ini & data crap****
print ''
name = 'WASP-53'        #default system
e1 = 0.01                #default e
if len(sys.argv) >= 2:
    name = str(sys.argv[1])
if len(sys.argv) >= 3:
    e1 = float(sys.argv[2])

#            name      M    R     a1     m1   r1    a2   e2   m2
data_bank=[('WASP-53',0.85,0.81,0.04106,2.13,1.074,3.39,0.829,16),
           ('WASP-81',1.07,1.28,0.03908,0.725,1.422,2.441,0.5667,57.3),
           ('HAT-P-13',1.22,1.56,0.04275,0.851,1.28,1.189,0.691,15.2)] #0.4383

for i in xrange(0,len(data_bank)):
    if name == data_bank[i][0]:
        params = data_bank[i]
        break

#********************
#conversion factors & constants
#J2S_M = 0.0009543       #mass Jupiter -> mass sun
J2S_R = 0.10045         #radius Jupiter -> radius sun
#J2AU_R = 0.0004673195   #radius Jupiter->AU
#c=10065.2               #speed of light AU/(yr/2pi)
#G=1                     #AU^3/Ms/(yr/2pi)^2

#SI
G=6.67*10**(-11)        #G in SI m^3/kg/s^2
J2kg_M = 1.898*10**27   #jupiter mass -> kg
S2kg_M = 1.989*10**30   #solar mass ->kg
J2m_R = 69911000        #jupiter radius -> meters
AU2m_a = 149597871000   #AU -> meters
c=299792458             #c in m/s

#star
M = params[1]*S2kg_M
#R = params[2]

#inner planet
a1 = params[3]*AU2m_a       #convert to meters from AU
m1 = params[4]*J2kg_M       #convert to kg from Jupiter mass
r1 = params[5]*J2m_R        #convert to meters from Jupiter radius

#outer planet
a2 = params[6]*AU2m_a       #convert to meters from AU
e2 = params[7]
m2 = params[8]*J2kg_M       #convert to kg from Jupiter mass

#choosing Q_inner
rp = params[5]*J2S_R        #For Q only - convert to solar radius from Jupiter radius.
if rp > 2*0.009156 and rp < 0.1:
    Q1 = 2.2e4
elif rp >= 0.1:
    Q1 = 5.4e4
else:
    Q1 = 40;

#precalculate terms that you can
#secular precession terms
ec2_inv = 1.0/(1.0 - e2*e2)
alpha = a1/a2
alpha3 = alpha**3
ab_3 = a1**3
nb = np.sqrt(G*(M+m1)/ab_3)
nb_3 = nb**3
nc = np.sqrt(G*(M+m2)/(a2**3))
cosdw = 1   #could also be -1
term1 = (1 + 4*e2*e2)*ec2_inv
#tidal terms
a5R5 = (a1/r1)**5
a_c2 = (a1/c)**2
e1_2 = e1*e1

#calculate k2 for each e1
num_points = 100
eb = np.zeros(num_points)
k2b = np.zeros(num_points)
for i in xrange(0,num_points):
    e1 = 0.04*i/num_points+0.001
    eb2_inv = 1.0/(1.0 - e1*e1)
    eb2_inv5 = eb2_inv**(-5)
    f2 = eb2_inv5*(1 + 1.5*e1_2 + 0.125*e1_2*e1_2)
    wb_GR = 3*nb_3*eb2_inv*a_c2
    wb_s = 0.75*nb*(m2/M)*alpha3*ec2_inv**(1.5)*(1 - 1.25*alpha*(e2/e1)*ec2_inv*cosdw)
    wc_s = 0.75*nc*(m1/M)*alpha*alpha*ec2_inv*ec2_inv*(1 - 1.25*alpha*(e1/e2)*term1*cosdw)
    numerator = 2*(wc_s - wb_s - wb_GR)*a5R5
    denominator = 15*(M/m1)*f2*nb + nb_3*ab_3*eb2_inv*eb2_inv/(G*m1)
    k2b[i] = numerator/denominator
    eb[i] = e1

plt.plot(k2b,eb, label='Calculated values')
plt.xlim([0,1.0])
plt.xlabel('k$_{2,b}$', fontsize=15)
plt.ylabel('e$_b$', fontsize=15)
plt.title(name)

#eq.16 - fit
if name == 'HAT-P-13':
    k2_fit = np.zeros(num_points)
    e_fit = np.zeros(num_points)
    for i in xrange(0,num_points):
        k2 = 1.5*i/num_points
        e_fit[i] = 0.0334 - 0.0985*k2 + 0.188*k2*k2 - 0.184*k2*k2*k2 + 0.069*k2*k2*k2*k2
        k2_fit[i] = k2
    plt.plot(k2_fit,e_fit, 'r',label='Batygin fit')
    plt.legend(loc='upper left',prop={'size':10})

plt.show()