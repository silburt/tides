#Purpose of this macro is to plot each system on the y-axis and a quantity on the x-axis,
#with probably a color gradient for a third variable. E.g. plot time on x-axis, and eccentricity
#is the color scheme (deep red = large eccentricity, light red = low eccentricity).

import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import matplotlib.cm as cm

files = glob.glob("globalplot/*.txt")
N_files = len(files)  #number of files we're dealing with
yaxisval=np.zeros(N_files)

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=int(sys.argv[1])   #x-axis variable
arg2=int(sys.argv[2])   #color coding variable

ranger=6
if arg2 == 1:
    xmax=100000000.

#for i in xrange(0,N_files):
for i in xrange(0,ranger):
    fos = open(''+files[i], 'r')
    header = fos.readline()
    header = header.split(",")
    files[i]=header[0]
    Ms = float(header[1])
    Rs = float(header[2])
    N = int(header[3])
    tide_delay = float(header[4])
    rp = np.zeros(N)
    mp = np.zeros(N)
    P = np.zeros(N)
    Qp = np.zeros(N)
    mig = np.zeros(N)
    for m in range(0,N):
        header = fos.readline()
    if arg2 == 1:
        data = np.loadtxt(fos, delimiter="	")
        for k in xrange(1,N):
            p=data[k::N]
            for j in xrange(0,k):
                q=data[j::N]
                if(abs(p[-1,3]/(2*q[-1,3]) - 1) < 0.05):
                    val = abs(2.0 - p[:,3]/q[:,3])       #Proximity from 2:1 res
                    gradient = 0 + (val[100:-1] - 0.03)*(1-0)/(-0.03)   #normalize between 0 and 1
                    yval=np.zeros(len(gradient)) + 0.2*i + 0.2
                    yaxisval[i] = 0.2*i + 0.2
                    plt.scatter(p[100:-1,arg1], yval, c=gradient, cmap=cm.coolwarm, lw=0)

#print yaxisval, files
plt.yticks(yaxisval[0:ranger], files[0:ranger], rotation='horizontal')
#fig, ay = plt.subplots()
#labels = [item.get_text() for item in ay.get_yticklabels()]
#labels = files
#ay.set_yticklabels(labels)

#plt.ylim([0.,0.11])
#plt.ylim([0.2025,0.2075])
plt.xlim([0,xmax])
#plt.title(''+name)
#plt.xlabel('' + names[arg1])
#plt.ylabel('' + names[arg2])
#plt.legend(loc='upper left',prop={'size':10})
plt.show()