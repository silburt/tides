import sys
import matplotlib.pyplot as plt
import numpy as np
import math

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
file_name='Kepler-85good.txt'

#Get basic system params from header of file
fos = open('xyplotter/'+file_name, 'r')
data = np.loadtxt(fos, delimiter="	")

N=4         #Number of planets
a=data[0::N]
b=data[1::N]
c=data[2::N]
d=data[3::N]

nchunks = 5
colors = ['black','brown','green','yellow','red']

p = a
steps = p.shape[0]/nchunks
for i in xrange(0,nchunks):
    plt.plot(d[i*steps:(i+1)*steps-1,1],d[i*steps:(i+1)*steps-1:,2], color=''+colors[i], alpha = 1.3 - 0.1*i , label='time chunk'+str(i+1))
    plt.plot(c[i*steps:(i+1)*steps-1,1],c[i*steps:(i+1)*steps-1:,2], color=''+colors[i], alpha = 1.3 - 0.1*i)
    plt.plot(b[i*steps:(i+1)*steps-1,1],b[i*steps:(i+1)*steps-1:,2], color=''+colors[i], alpha = 1.3 - 0.1*i)
    plt.plot(a[i*steps:(i+1)*steps-1,1],a[i*steps:(i+1)*steps-1:,2], color=''+colors[i], alpha = 1.3 - 0.1*i)

plt.xlabel('x (AU)')
plt.ylabel('y (AU)')
plt.xlim([-0.3,0.3])
plt.ylim([-0.3,0.3])
plt.legend(loc='upper left',prop={'size':10})
plt.title(''+file_name)
plt.show()