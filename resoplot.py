#The purpose of this macro is to plot the final period ratios of the Kepler systems and compare it to the original values.
import sys
import numpy as np
import matplotlib.pyplot as plt

arg1='2'
arg2='1'
thresh=0.05
simpath='saved_runs/round5_Feb27'

systems=np.genfromtxt('reso/'+arg1+':'+arg2+'_systems.txt', dtype=("|S10")) #system names
N_sys = len(systems)
K85i = 0
K223i = 0
for i in xrange(0,N_sys):
    if systems[i] == 'Kepler-223':
        K223i = i
    if systems[i] == 'Kepler-85':
        K85i = i
systems = np.delete(systems,[K223i,K85i])    #Remove these systems for now
N_sys -= 2

data = np.genfromtxt('planets.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float), filling_values="0", skip_header=34)
n_data = data.shape[0]
PRobs=np.zeros(0)       #obs period ratios
PRsim=np.zeros(0)       #sim period ratios
e_in = np.zeros(0)      #eccentricity of inner planet < 0.06
e_inL = np.zeros(0)     #eccentricity of inner planet > 0.06
time=np.zeros(N_sys)
Navg = 100
Pavg=np.zeros(Navg-1)

for k in xrange(0,N_sys): #find period ratio of planets
    name = systems[k]
    print name
    found = 0
    index = 0
    while found == 0:    #look at obs planets
        if name == data[index][1]:
            N = data[index][4]
            found = 1
        else:
            index +=1
    for i in xrange(index, index+N-1):
        for j in xrange(i+1, index+N):
            ratio = float(data[j][5])/float(data[i][5])
            arg = ratio - float(arg1)/float(arg2)
            if abs(arg) < thresh and arg > 0.: #Is each i,j in *outer* res
                PRobs = np.append(PRobs,ratio)
                outer = j-index
                inner = i-index
    fos = open(simpath+'/orbits_'+name+'.txt','r')
    lines = fos.readlines()
    length = len(lines)
    for l in xrange(1,Navg):
        temp = lines[length-l*N+outer]  #outer planet avg period value over last Navg orbits
        tempP = temp.split("\t")
        Pavg[l-1] = float(tempP[3])
    outp = np.mean(Pavg)
    for l in xrange(1,Navg):
        temp = lines[length-l*N+inner]
        tempP = temp.split("\t")
        Pavg[l-1] = float(tempP[3])     #inner planet avg period value over last Navg orbits
    inp = np.mean(Pavg)
    simratio = outp/inp
    PRsim = np.append(PRsim,simratio)
    tempe = lines[length-N+inner]
    tempe = tempe.split("\t")
    time[k] = float(tempe[0])/1000000.
    e = float(tempe[2])
    if e > 0.05:
        e = 0.05001
        e_inL = np.append(e_inL,e)
    else:
        e_in = np.append(e_in,e)   #eccentricity of inner planet (estimate of tides)

fos.close()
binwidth=0.005
cdf=1
#cdf
if cdf == 1:
    plt.hist(PRobs - 2.0, color='black', linewidth=2, bins=np.arange(0., 0.06 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label='Observed $\Delta$ from 2:1 MMR')
    plt.hist(PRsim - 2.0, color='green', linewidth=2, bins=np.arange(0., 0.06 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label='Simulated $\Delta $ from 2:1 MMR after '+str(round(max(time)))+'Myr')
    plt.hist(PRobs - PRsim, color='yellow', linewidth=2, bins=np.arange(0., 0.06 + binwidth, binwidth), histtype='step',cumulative='true', normed='true',label='$\Delta_{Obs}$ - $\Delta_{Sim}$ (system-by-system basis)')
    #plt.hist(e_in, alpha=0.5, color='orange', bins=np.arange(0., 0.06 + binwidth, binwidth), normed='true', label='Simulated e$_{inner}$')
    #plt.hist(e_inL, alpha=0.5, color='yellow', bins=np.arange(0., 0.06 + binwidth, binwidth), label='Simulated e$_{inner}$ > 0.05')
else:
    plt.hist(PRobs - 2.0, color='black', bins=np.arange(0., 0.06 + binwidth, binwidth), label='Observed $\Delta$ from 2:1 MMR')
    plt.hist(PRsim - 2.0, color='green', alpha=0.8, bins=np.arange(0., 0.06 + binwidth, binwidth), label='Simulated $\Delta$ from 2:1 MMR after '+str(round(max(time)))+'Myr')
    plt.hist(PRobs - PRsim, color='yellow', alpha=0.4, bins=np.arange(0., 0.06 + binwidth, binwidth), label='$\Delta_{Obs}$ - $\Delta_{Sim}$ (system-by-system basis)')

plt.xlabel('P$_2$/P$_1$')
plt.ylabel('counts, total='+str(len(PRsim)))
plt.title('Period ratios of Kepler planets close to 2:1 MMR')
plt.legend(loc='upper left',prop={'size':10})
plt.show()




