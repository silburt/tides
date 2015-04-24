import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import math
import matplotlib.cm as cm
pi = math.pi

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
file_name=str(sys.argv[1])
arg1=int(sys.argv[2])
arg2=int(sys.argv[3])
arg3 = -1
arg4 = 0
if len(sys.argv) == 5:
    arg3 = int(sys.argv[4])
elif len(sys.argv) == 6:
    arg3 = int(sys.argv[4])
    arg4 = int(sys.argv[5])


#Load numerical data
#names=['time (years)','Semi-Major Axis (AU)','Eccentricity','Period (Days)','arg. of peri','Mean Anomaly','Eccentric Anomaly','Mean Longitude (lambda)','Resonant Angle (phi = 2*X2 - X1 - w1)','Resonant Angle2 (phi2 = 2*X2 - X1 - w2)','Libration Timescale (order of mag.)','Period Ratio (P$_{i+1}$/P$_{i}$) - j/(j+1)','Resonance Plot','G/G0 - 1']
names=['time (years)','r','vr','tidetau_e','tidetau_a','dvx','dvy','-dvx/(2.*tidetau_a[i])','-dvy/(2.*tidetau_a[i])','-2./tidetau_e[i]*vr*dx/r','-2./tidetau_e[i]*vr*dy/r','-dvx/(2.*tidetau_a[i]) - 2./tidetau_e[i]*vr*dx/r','-dvy/(2.*tidetau_a[i]) - 2./tidetau_e[i]*vr*dy/r','dxz','dz']
colors=['b','g','m','r','c','y']
fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=",")
N=2

for i in range(0,N): #range(0,N) only goes to N-1
    p=data[i::N]
    plt.plot(p[arg4:arg3,arg1], p[arg4:arg3,arg2], 'o'+colors[i], label='planet '+str(i+1), markeredgecolor='none')

plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.legend(loc='upper left',prop={'size':10})
plt.show()
