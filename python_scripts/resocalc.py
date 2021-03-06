#Purpose of this macro is to find the systems within X percent of the specified MMR

import sys
import numpy as np
import matplotlib.pyplot as plt

arg1=int(sys.argv[1])   #(j+1) resonance
arg2=int(sys.argv[2])   #j resonance

#File for writing
output = open('../reso/'+str(arg1)+':'+str(arg2)+'_systems.txt','w')
outputfull = open('../reso/full/'+str(arg1)+':'+str(arg2)+'_systems_fulldetail.txt','w')

data=np.genfromtxt('../planets.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float), filling_values="0", skip_header=34)

Pratio=np.zeros(1)

nlines=data.shape[0]
count=0
Ndsys=0
thresh=0.06
o_index=0
large_mass=0
skip=0
while (count <= nlines - 1):
    name = data[count][1]   #starting place of given system
    Ndsys = count
    while (data[Ndsys][1] == name):
        Ndsys += 1           #eNding place of given system
        if Ndsys > nlines - 1:
            break
    #for i in xrange(count,Ndsys-1):
        #if data[i][22] > 0.2 or data[i][19] > 100. :   #comment out small mass requirement??
            #large_mass = 1
    if not 'TEST' in name and large_mass != 1:  #look for resonances within system if small masses/radii & not a test case
        for i in xrange(count, Ndsys-1):
            for j in xrange(i+1, Ndsys):
                ratio = float(data[j][5])/float(data[i][5])
                arg = ratio - float(arg1)/float(arg2)
                if abs(arg) < thresh and arg > 0.: #Is each i,j in *outer* res
                    Pratio= np.vstack([Pratio,ratio])
                    length = len(data[0])                         #Full output for detail
                    for k in xrange(0,length - 2):
                        outputfull.write(str(data[i][k])+',')
                    outputfull.write(str(data[i][length-1])+','+str(i-count)+'\n')
                    for k in xrange(0,length - 2):
                        outputfull.write(str(data[j][k])+',')
                    outputfull.write(str(data[j][length-1])+','+str(j-count)+'\n')
                    if skip != 1:
                        output.write(data[i][1]+'\n')
                        skip = 1
    count = Ndsys
    large_mass = 0
    skip=0

output.close()
Pratio = np.delete(Pratio, 0,0) #delete first row that was created earlier
#cdf = np.cumsum(Pratio)
binwidth=0.01
center = float(arg1)/float(arg2)
plt.hist(Pratio, bins=np.arange(center - thresh, center + thresh + binwidth, binwidth), color='gray')
plt.plot([2.06, 2.06], [0,10], 'r--', linewidth = 2)
plt.xlabel('P$_2$/P$_1$')
plt.ylabel('counts')#, total='+str(len(Pratio)))
#plt.title('Histogram of Planets close to Resonance')
plt.xlim([center - thresh, center + thresh])
plt.show()




