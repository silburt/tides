#The purpose of this macro is to calculate the (instantaneous) equilibrium eccentricity of a planet, which is eq. 36 from Mardling (2007).

import sys
import numpy as np
import matplotlib.pyplot as plt
import math

#ini & data crap****
print ''
name = 'WASP-53'        #default system
if len(sys.argv) >= 2:
    name = str(sys.argv[1])

data_bank=[('WASP-53',0.85,0.81,0.04106,2.13,1.074,3.39,0.829,16),
           ('WASP-81',1.07,1.28,0.03908,0.725,1.422,2.441,0.5667,57.3)]

for i in xrange(0,2):
    if name == data_bank[i][0]:
        params = data_bank[i]
        break

save_output = 1
#********************
#star
M = params[1]
R = params[2]

#conversion factors
J2S_M = 0.0009543       #mass Jupiter -> mass sun
J2S_R = 0.10045         #radius Jupiter -> radius sun
J2AU_R = 0.0004673195   #radius Jupiter->AU

#inner planet
a1 = params[3]
m1 = params[4]*J2S_M
r1 = params[5]*J2AU_R

#choosing Q_inner
rp = params[5]*J2S_R
if rp > 2*0.009156 and rp < 0.1:
    Q1 = 2.2e4
elif rp >= 0.1:
    Q1 = 5.4e4
else:
    Q1 = 40;

if len(sys.argv) >=3:
    k1 = float(sys.argv[2])

#outer planet
a2 = params[6]
e2 = params[7]
m2 = params[8]*J2S_M

#constants
c=10065.2 #speed of light AU/(yr/2pi)
G=1

#e_eq (Eq. 36 from Mardling, 2007)
ep_c = np.sqrt(1.0 - e2*e2)
n1_sq = (G*(M + m1))/(a1*a1*a1)
gamma = 4*(a1*a1*n1_sq/(c*c))*(M/m2)*((a2/a1)**3)
num = (5.0/4.0)*(a1/a2)*e2/(ep_c*ep_c)
den = abs(1.0 - np.sqrt(a1/a2)*(m1/m2)/ep_c + gamma*(ep_c**3))
e_eq = num/den
print 'e_eq for WASP-53b is',e_eq

#solve for k2 of inner planet (Eq. 16 from Batygin 2009)
max = 1000
convergence = 0
min_result = 0.334
k2 = 1      #defualt k value for inner planet
for i in xrange(0,max):
    k2 = i*1.5/max
    result = 0.0334 - 0.0985*k2 + 0.188*(k2**2) - 0.184*(k2**3) + 0.069*(k2**4)
    if abs(result - e_eq) < 0.001:
        print '(Batygin) Calculated k2 of inner planet is', k2
        convergence = 1
        break
if convergence == 0:
    print '(Batygin) No convergence for k2.'

#tidal timescale of inner planet
R5a5 = (r1/a1)**5
GM3a3 = (M/a1)**1.5
edot_e = (21/2.)*(k2/Q1)*GM3a3*R5a5/m1
tau_e = 1.0/edot_e/1e6
print 'tau_e for WASP-53b is',tau_e,' Myr.'

if save_output == 1:
    text_file = open("eqcalc/eqcalc.txt", "a")
    text_file.write("%s:\t%f\t%f\t%.1f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (name,m1,r1,Q1,m2,a2,e2,e_eq,tau_e,k2))
    text_file.close()
