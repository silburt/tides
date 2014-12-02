import sys
import matplotlib.pyplot as plt
import numpy as np
import math

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=int(sys.argv[1])
arg2=int(sys.argv[2])
names=['time (years)','Semi-Major Axis (AU)','Eccentricity','inclination','Long. Asc. Node','arg. of peri','Mean Longitude','Period (Days)','mean anomaly','period ratio']

fos = open('runs/orbits_notides_Kepler92.txt', 'r')
data = np.loadtxt(fos, delimiter="	")
p1 = data[0::2] #2 planets, start at 0th entry, get every 2nd entry
p2 = data[1::2]

if arg2 == 9:
    plt.plot(p1[:,arg1], p2[:,7]/p1[:,7], 'b')
else:
    plt.plot(p1[0:2000,arg1], p1[0:2000,arg2], 'b')
    plt.plot(p2[0:2000,arg1], p2[0:2000,arg2], 'g')

plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.show()