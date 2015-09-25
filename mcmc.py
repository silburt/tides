#Purpose of this macro is to find the optimum parameters (via MCMC) for a given system which maximize the span*k2_at_half_max value. Input known parameters of system first.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import random

#conversion factors & constants - G=1 units
J2S_M = 0.0009543       #mass Jupiter -> mass sun
J2S_R = 0.10045         #radius Jupiter -> radius sun
J2AU_R = 0.0004673195   #radius Jupiter->AU
c=10065.2               #speed of light AU/(yr/2pi)
G=1                     #AU^3/Ms/(yr/2pi)^2

def conversion_factors(data, sigma_step_factor,upper_bound,lower_bound):
    sigma_step = data/sigma_step_factor
    data[2] *= J2S_M
    data[3] *= J2AU_R
    data[6] *= J2S_M
    sigma_step[2] *= J2S_M
    sigma_step[3] *= J2AU_R
    sigma_step[6] *= J2S_M
    upper_bound[2] *= J2S_M
    upper_bound[3] *= J2AU_R
    upper_bound[6] *= J2S_M
    lower_bound[2] *= J2S_M
    lower_bound[3] *= J2AU_R
    lower_bound[6] *= J2S_M
    return data, sigma_step, upper_bound, lower_bound

#*************************FUNCTIONS************************************
def calc_span_k2(data,param_test,N_params,param_index,num_points_e):
#ini params
    eb = np.zeros(0)     #reset values
    k2b = np.zeros(0)
    k2prev = 0
    minlength = 30
    max_e = 0.15
#calcs
    for i in xrange(0,N_params):
        data[param_index[i]] = param_test[i]
    M,a1,m1,r1,a2,e2,m2 = data
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
        elif (k2prev < 0 or k2prev > 1.5) and len(k2b) > minlength:
            #print 'Just calcing 0 > k2 > 1.5 values. Breaking early at iteration',i
            break
        k2prev = k2val
    if len(k2b) > minlength:  #make sure there's actually enough useable values
        maximum = max(eb)
        minimum = min(eb)
        half = 0.5*(maximum - minimum) + minimum
        for k in xrange(0,num_points_e):
            if half <= eb[k]:
                return (maximum - minimum)*k2b[k]
        print '!ERROR! Couldnt find half max! Need to debug. Exiting.'
        exit()
    else:
        print 'Not enough useable values to calculate span. Span_k2 = 0'
        return 0

#***************************SETUP**************************************
#****DATA************
#ini values  M=0   a1=1   m1=2 r1=3 a2=4  e2=5  m2=6   name=7
data_bank=[(1.22,0.04275,0.851,1.28,1.189,0.691,15.2,'HAT-P-13'),
           (1.126,0.080,0.442,1.18,2.39,0.21,3.71,'KELT-6')]
sigma_step_factor=100      #how big is the jump size relative to the data value?

system = 'HAT-P-13'

data = np.zeros(7)
for i in xrange(0,len(data_bank)):
    if system == data_bank[i][7]:
        for j in xrange(0,7):
            data[j] = data_bank[i][j]
        break

#upper/lower bounds on data
upper_bound=[10,data[4],20,2,data[4]*5,0.99,20]
lower_bound=[0.1,0.01,0,0.01,data[1],0,0]

data, sigma_step,upper_bound,lower_bound = conversion_factors(data, sigma_step_factor,upper_bound,lower_bound)

#****INI PARAMS*****
N_params = 4
param_index = [0,2,3,6]     #index of variables that are going to be varied in the MCMC. M=0, a1=1, etc.
num_curve_points = 500      #resolution of curve with which the span and k2_at_half_max is calc'd from
N_iterations = 1000
n_jumps = 0

#ini arrays
param_results = np.zeros((N_iterations,N_params))     #store the results of each parameter
span_k2 = np.zeros(N_iterations)
param_test = np.zeros(N_params)
for i in xrange(0,N_params):
    param_test[i] = data[param_index[i]]
    param_results[0,i] = param_test[i]

#random number generator
random.seed(1)               #if no arg. it uses the current time as seed

#***************************MCMC**************************************
#ini MCMC
span_k2_max = calc_span_k2(data,param_test,N_params,param_index,num_curve_points)
for i in xrange(1,N_iterations):
    for j in xrange(0,N_params):    #generate new random parameters
        rnd = random.gauss(0,1)
        #print 'vals',param_test[j],rnd,sigma_step[j],rnd*sigma_step[j],param_results[i-1,j]
        param_test[j] = param_results[i-1,j] + rnd*sigma_step[j]
        while param_test[j] > upper_bound[param_index[j]] or param_test[j] < lower_bound[param_index[j]]:
            param_test[j] = param_results[i-1,j] + random.gauss(0,1)*sigma_step[j]
    span_k2_current = calc_span_k2(data,param_test,N_params,param_index,num_curve_points)
    if span_k2_current > span_k2_max:
        span_k2_max = span_k2_current
        param_results[i,:] = param_test
        n_jumps += 1
    else:
        p = np.exp(-50*(span_k2_max - span_k2_current))
        rand = random.random()
        #print 'jump?',p, rand, span_k2_max, span_k2_current
        if rand > p:
            span_k2_max = span_k2_current
            param_results[i,:] = param_test
            n_jumps += 1
        else:
            param_results[i,:] = param_results[i-1,:]
    span_k2[i] = span_k2_current    #keep track of the evolution of span_k2
    if float(i)/100. - i/100 == 0:
        print 'iteration', i

