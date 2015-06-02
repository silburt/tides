#Purpose of this macro is to display the distribution of at what eccentricity do systems go unstable at. I.e. Kepler_e_coll.txt has printed values where, given an initial eccentricity (to each planet) of X, does the system go unstable in **2 million years** (time is important)?

import sys
import numpy as np
import matplotlib.pyplot as plt

N = 5
e_count = np.zeros(N)   #counting e_i=0.2,0.4,0.6,0.75
x = np.arange(N)

def options(x):
    return {
        0.2 : 0,
        0.3 : 1,
        0.4 : 2,
        0.6 : 3,
        0.75 : 4,
    }[x]

data=np.genfromtxt('emax_dyn.txt', delimiter=',', dtype=("|S10",float,float))
nlines=data.shape[0]
N_sys = 34 #real total = 34, but some systems I'm excluding

for i in xrange(0,nlines):
    num = float(data[i][1])
    e_count[options(num)] += 1

x_arr = [0,1,2,3,4]
fig, ax = plt.subplots()
ax.plot(x_arr, e_count/N_sys, '.', ms=15, color='black')
plt.errorbar(x_arr,e_count/N_sys,yerr=np.sqrt(e_count)/N_sys, color='black',linestyle='none', linewidth=2)
ax.set_xlim([-0.2, 4.2])
ax.set_ylim([0,1])
handles, labels = ax.get_legend_handles_labels()
ax.set_xticks(x)
ax.set_xticklabels( ('0.2', '0.3', '0.4', '0.6', '0.75') )
ax.set_xlabel('$e_i$ of each planet in system', fontsize=16)
ax.set_ylabel('Fraction of Systems With Collision in 2Myr', fontsize=16)
plt.show()



