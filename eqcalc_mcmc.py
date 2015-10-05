#Purpose of this macro is to find the optimum parameters (via MCMC) for a given system which maximize the span*k2_at_half_max value. Input known parameters of system first.

import sys
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
import random

#conversion factors & constants - G=1 units
J2S_M = 0.0009543       #mass Jupiter -> mass sun
J2S_R = 0.10045         #radius Jupiter -> radius sun
J2AU_R = 0.0004673195   #radius Jupiter->AU
c=10065.2               #speed of light AU/(yr/2pi)
G=1                     #AU^3/Ms/(yr/2pi)^2

def conversion_factors(data,sigma_factor,upper_bound,lower_bound,param_index,N_params):
    sigma_step = np.zeros(N_params)
    for i in xrange(0,N_params):
        PI = param_index[i]
        sigma_step[i] = data[PI]/sigma_factor[PI]
        if PI == 2 or PI == 6:  #convert to G=1 units
            sigma_step[i] *= J2S_M
        elif PI == 3:
            sigma_step[i] *= J2AU_R
    data[2] *= J2S_M
    data[3] *= J2AU_R
    data[6] *= J2S_M
    upper_bound[2] *= J2S_M
    upper_bound[3] *= J2AU_R
    upper_bound[6] *= J2S_M
    lower_bound[2] *= J2S_M
    lower_bound[3] *= J2AU_R
    lower_bound[6] *= J2S_M
    return data, sigma_step, upper_bound, lower_bound

def conversion_back_to_norm(params,param_index,N_params):
    for i in xrange(0,N_params):
        if param_index[i] == 2 or param_index[i] == 6:
            params[:,i] /= J2S_M
        elif param_index[i] == 3:
            params[:,i] /= J2AU_R
    return params

#*************************FUNCTIONS************************************
def calc_span_k2(data,param_test,N_params,param_index):
#ini params
    eb = np.zeros(0)     #reset values
    k2b = np.zeros(0)
    k2prev = 0
    span_minvals = 30
    num_curve_points = 1000      #resolution of curve with which the span and k2_at_half_max is calc'd from
    max_e = 0.2
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
    for i in xrange(0,num_curve_points):
        e1 = max_e*i/num_curve_points+1e-5
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
        elif (k2prev < 0 or k2prev > 1.5) and len(k2b) > span_minvals:
            #print 'Just calcing 0 > k2 > 1.5 values. Breaking early at iteration',i
            break
        k2prev = k2val
    if len(k2b) > span_minvals:  #make sure there's actually enough useable values
        maximum = max(eb)
        minimum = min(eb)
        half = 0.5*(maximum - minimum) + minimum
        for k in xrange(0,num_curve_points):
            if half <= eb[k]:
                return (maximum - minimum)*k2b[k]
        print '!ERROR! Couldnt find half max! Need to debug. Exiting.'
        exit()
    else:
        print 'Not enough useable values to calculate span. Skip Calculation.'
        return -1

#***************************SETUP**************************************
#****INI PARAMS*****
system = 'Neptune'
N_iterations = 10000
N_params = 3
param_index = [0,1,6]     #index of variables that are going to be varied in the MCMC. M=0, a1=1, etc.

#****DATA************
#ini values  M=0   a1=1   m1=2 r1=3 a2=4  e2=5  m2=6   name=7
data_bank=[(1.22,0.04275,0.851,1.28,1.189,0.691,15.2,'HAT-P-13'),
           (1.0,0.1,0.05,0.352,1.0,0.7,2.0,'Neptune'),
           (1.126,0.080,0.442,1.18,2.39,0.21,3.71,'KELT-6')]
sigma_factor=[50,50,50,50,50,50,30]    #how big is the jump size relative to the data value?

data = np.zeros(7)
for i in xrange(0,len(data_bank)):
    if system == data_bank[i][7]:
        for j in xrange(0,7):
            data[j] = data_bank[i][j]
        break
data_copy = data

#upper/lower bounds on data
#M=0   a1=1   m1=2 r1=3 a2=4  e2=5  m2=6
upper_bound=[10,data[4],20,2,data[4]*5,0.99,20]
lower_bound=[0.01,0.01,0,0.01,data[1],0,0]

data,sigma_step,upper_bound,lower_bound = conversion_factors(data,sigma_factor,upper_bound,lower_bound,param_index,N_params)

#ini arrays
param_results = np.zeros((N_iterations,N_params))     #store the results of each parameter
span_k2 = np.zeros(N_iterations)

param_test = np.zeros(N_params)
for i in xrange(0,N_params):
    param_test[i] = data[param_index[i]]
    param_results[0,i] = param_test[i]

#final ini params
random.seed(2)               #if no arg. it uses the current time as seed
n_jumps = 0
prob_weight = -2
print_increment = 500

#***************************MCMC**************************************
#ini MCMC
span_k2_max = calc_span_k2(data,param_test,N_params,param_index)
for i in xrange(1,N_iterations):
    for j in xrange(0,N_params):    #generate new random parameters
        rnd = random.gauss(0,1)
        param_test[j] = param_results[i-1,j] + rnd*sigma_step[j]
        while param_test[j] > upper_bound[param_index[j]] or param_test[j] < lower_bound[param_index[j]]:
            param_test[j] = param_results[i-1,j] + random.gauss(0,1)*sigma_step[j]
    span_k2_current = calc_span_k2(data,param_test,N_params,param_index)
    if span_k2_current < 0:
        span_k2[i] = span_k2[i-1]
        param_results[i,:] = param_results[i-1,:]
    else:
        if span_k2_current > span_k2_max:
            span_k2_max = span_k2_current
            param_results[i,:] = param_test
            n_jumps += 1
        else:
            #erf = (span_k2_max - span_k2_current)/span_k2_max
            erf = span_k2_max - span_k2_current
            prob_weight = -50
            p = np.exp(prob_weight*erf)
            rand = random.random()
            #print 'jump?',p, rand, span_k2_max, span_k2_current
            if rand > p:
                span_k2_max = span_k2_current
                param_results[i,:] = param_test
                n_jumps += 1
            else:
                param_results[i,:] = param_results[i-1,:]
        span_k2[i] = span_k2_current    #keep track of the evolution of span_k2
    if float(i)/float(print_increment) - i/print_increment == 0:
        print 'iteration', i

param_results = conversion_back_to_norm(param_results,param_index,N_params)  #convert back to orig
names = ['M','a1','m1','r1','a2','e2','m2']
units = ['M$_{sun}$','AU','M$_{Jup}$','R$_{Jup}$','AU','','M$_{Jup}$']

#plotting and results
fig, a = plt.subplots(nrows=4, ncols=1, figsize=(10,10))
fig.subplots_adjust(left=0.12, right=0.94)
title = system+': MCMC Results, varying '+str(N_params)+' params'
fig.text(0.5, 0.95, title, ha='center', va='center', rotation='horizontal', fontsize=20)
print
print 'Best Fit Parameters for System '+system
burnin = int(N_iterations*0.1)      #first 10% = burn in?
length = N_iterations - burnin
for i in xrange(0,N_params):
    a[i].plot(param_results[:,i])
    a[i].set_ylabel(names[param_index[i]]+' ('+units[param_index[i]]+')', fontsize=15)
    a[i].ticklabel_format(useOffset=False)
    result = param_results[burnin:N_iterations-1,i]
    result.sort()
    med = result[0.5*length]
    high = result[0.84*length] - med
    low = result[0.16*length]
    a[i].plot([0,N_iterations],[med,med],'--',color='black')
    print names[param_index[i]],'='+str(med)+' + '+str(high)+' - '+str(low)+' '+units[param_index[i]]
a[N_params].plot(span_k2)
a[N_params].set_xlabel('Iteration', fontsize=15)
a[N_params].set_ylabel('span*k2', fontsize=15)
print 'max span*k2 = ',span_k2_max
print '%_jumps =',n_jumps/float(N_iterations)

fig.savefig('eqcalc/eqmcmc/'+system+'_mcmc_N='+str(N_iterations)+'.png')
plt.show()
