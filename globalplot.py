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
#arg1=int(sys.argv[1])   #x-axis variable
#arg2=int(sys.argv[2])   #color coding variable
arg1=0
arg2=1

if arg2 == 1:
    xmax=100000000.

#for i in xrange(0,N_files):
for i in xrange(0,N_files):
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
                if(abs(p[-1,3]/(2*q[-1,3]) - 1) < 0.05): #check last elements
                    out = p[10::10]    #start at 10th, and take every 10th output for brevity
                    length=len(out[:,8])
                    x = np.zeros(length)
                    y = np.zeros(length)
                    for m in xrange(0,length):
                        x[m] = out[m,2]*math.cos(out[m,8])
                        y[m] = out[m,2]*math.sin(out[m,8])
                    gradient = 0 + (out[:,8] - out[50,8])*(1-0)/out[50,8]
                    #val = abs(2.0 - p[:,3]/q[:,3])       #Proximity from 2:1 res
                    #gradient = 0 + (val[100:-1] - 0.03)*(1-0)/(-0.03)   #normalize between 0 and 1
                    yval=np.zeros(length) + 0.2*i + 0.2
                    yaxisval[i] = 0.2*i + 0.2
                    plt.scatter(out[:,arg1], yval, c=gradient, cmap=cm.coolwarm, lw=0)

#print yaxisval, files
plt.yticks(yaxisval[0:N_files], files[0:N_files], rotation='horizontal')
#fig, ay = plt.subplots()
#labels = [item.get_text() for item in ay.get_yticklabels()]
#labels = files
#ay.set_yticklabels(labels)

#plt.ylim([0.,0.11])
#plt.ylim([0.2025,0.2075])
plt.xlim([0,xmax])
cbar = plt.colorbar()
cbar.set_label('proximity from 2:1 resonance')
#plt.title(''+name)
#plt.xlabel('' + names[arg1])
#plt.ylabel('' + names[arg2])
#plt.legend(loc='upper left',prop={'size':10})
plt.show()