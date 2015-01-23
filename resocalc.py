import sys
import numpy as np
import matplotlib.pyplot as plt

arg1=int(sys.argv[1])   #(j+1) resonance
arg2=int(sys.argv[2])   #j resonance

#File for writing
output = open('reso/'+str(arg1)+':'+str(arg2)+'_systems.txt','w')
outputfull = open('reso/full/'+str(arg1)+':'+str(arg2)+'_systems_fulldetail.txt','w')

data=np.genfromtxt('planets.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float), filling_values="0", skip_header=34)

Pratio=np.zeros(1)

nlines=data.shape[0]
count=0
Ndsys=0
thresh=0.1
o_index=0
while (count <= nlines - 1):
    name = data[count][1]   #starting place of given system
    Ndsys = count
    while (data[Ndsys][1] == name):
        Ndsys += 1           #eNding place of given system
        if Ndsys > nlines - 1:
            break
    if not 'TEST' in name:  #look for resonances within system
        for i in xrange(count, Ndsys-1):
            for j in xrange(i+1, Ndsys):
                ratio = float(data[j][5])/float(data[i][5])
                if abs(ratio - float(arg1)/float(arg2)) < thresh: #Is each i,j in res
                    Pratio= np.vstack([Pratio,ratio])
                    length = len(data[0])                         #Full output for detail
                    for k in xrange(0,length - 2):
                        outputfull.write(str(data[i][k])+',')
                    outputfull.write(str(data[i][length-1])+'\n')
                    for k in xrange(0,length - 2):
                        outputfull.write(str(data[j][k])+',')
                    outputfull.write(str(data[j][length-1])+'\n\n')
                    output.write(data[i][1]+'\n')
    count = Ndsys

output.close()
Pratio = np.delete(Pratio, 0,0) #delete first row that was created earlier
#cdf = np.cumsum(Pratio)
binwidth=0.01
plt.hist(Pratio, bins=np.arange(round(100*min(Pratio))/100., round(100*max(Pratio))/100. + binwidth, binwidth))
plt.xlabel('P$_2$/P$_1$')
plt.ylabel('counts, total='+str(len(Pratio)))
plt.title('Histogram of Planets close to Resonance')
plt.show()




