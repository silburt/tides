#Purpose of this macro is to plot each system on the y-axis and a quantity on the x-axis,
#with probably a color gradient for a third variable. E.g. plot time on x-axis, and eccentricity
#is the color scheme (deep red = large eccentricity, light red = low eccentricity).

import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import glob
import matplotlib.cm as cm
pi = 3.14159265359

def res_index(mig,p,N):
    ini = np.amax(mig)
    fini = ini + 20000
    i_tide = 0
    f_tide = 0
    var = 0
    p=data[0::N]
    for i in xrange(0,len(p[:,0])):
        if p[i,0] > ini and var==0:
            i_tide = i
            var = 1
        if p[i,0] > fini:
            f_tide = i
            break
    if f_tide == i_tide:
        f_tide += 1
    return [i_tide,f_tide]

files = glob.glob("globalplot/*.txt")
N_files = len(files)  #number of files we're dealing with
yaxisval=np.zeros(N_files)

for i in xrange(0,N_files):
    fos = open(''+files[i], 'r')
    header = fos.readline()
    header = header.split(",")
    files[i]=header[0]
    Ms = float(header[1])
    Rs = float(header[2])
    N = int(header[3])
    tide_delay = float(header[4])
    mig = np.zeros(N)
    for m in range(0,N):
        header = fos.readline()
        header = header.split(",")
        mig[m] = float(header[5]) + float(header[5])/3. #damping time too
    data = np.loadtxt(fos, delimiter="	")
    for k in xrange(1,N):
        p=data[k::N]
        for j in xrange(0,k):
            q=data[j::N]
            ratio = p[-1,3]/(2*q[-1,3]) - 1
            if(abs(ratio) < 0.05) and ratio > -0.02:
                out = p[10::10]    #take every 10th output for brevity
                length=len(out[:,8])
                index = res_index(mig,p,N)
                i_min = index[0]
                i_max = index[1]
                for m in xrange(0,i_max+1):
                    if p[m,8] > pi:
                        p[m,8] -= 2.0*pi
                min = np.min(p[i_min:i_max,8])
                max = np.max(p[i_min:i_max,8])
                ref = (min + max)/2.
                for m in xrange(0,length):
                    if out[m,8] > ref + pi:
                        out[m,8] -= 2*pi
                    if out[m,8] < ref - pi:
                        out[m,8] += 2*pi
                gradient = abs(out[:,8] - ref)/pi
                yval=np.zeros(length) + 0.2*i + 0.2
                yaxisval[i] = 0.2*i + 0.2
                plt.scatter(out[:,0], yval, c=gradient, cmap=cm.coolwarm, vmin=0, vmax=1, lw=0)
                print 'completed object',i+1,' of ',N_files, ', System',files[i]

#print yaxisval, files

#fig, ay = plt.subplots()
#labels = [item.get_text() for item in ay.get_yticklabels()]
#labels = files
#ay.set_yticklabels(labels)

plt.yticks(yaxisval[0:N_files], files[0:N_files], rotation='horizontal')
xmax=500000000.
plt.xlim([0,xmax])
plt.ylim([0.15,0.2*N_files + 0.25])
cbar = plt.colorbar()
cbar.set_label('proximity from 2:1 resonance')
plt.title('Kepler planets in 2:1 resonance')
plt.show()