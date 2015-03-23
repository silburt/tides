#The purpose of this macro is to plot the final period ratios of the Kepler systems and compare it to the original values.
import sys
import numpy as np
import matplotlib.pyplot as plt

arg1='2'
arg2='1'
thresh=0.06
simpath='saved_runs/round8_Mar16Qpfac1'
ext = '_Qpfac1'

early = 1   #if = 1, plot just before tides start (ini conditions)

systems=np.genfromtxt('reso/'+arg1+':'+arg2+'_systems.txt', dtype=("|S10")) #system names
N_sys = len(systems)
N_remove = 4
Ki = np.zeros(N_remove)
for i in xrange(0,N_sys):
    if systems[i] == 'Kepler-223':
        Ki[0] = i
    if systems[i] == 'Kepler-85':
        Ki[1] = i
    if systems[i] == 'Kepler-31':
        Ki[2] = i
    if systems[i] == 'Kepler-11':
        Ki[3] = i
systems = np.delete(systems,Ki)    #Remove these systems for now
N_sys -= N_remove

data = np.genfromtxt('planets.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float), filling_values="0", skip_header=34)
n_data = data.shape[0]
PRobs=np.zeros(0)       #obs period ratios
PRsim=np.zeros(0)       #sim period ratios
e_in = np.zeros(0)      #eccentricity of inner planet < 0.06
e_inL = np.zeros(0)     #eccentricity of inner planet > 0.06
e0 = np.zeros(0)        #eccentricity of inner planet just before tides.
p0 = np.zeros(0)        #period ratio just before tides
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
    fos = open(simpath+'/orbits_'+name+''+ext+'.txt','r')
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
    #calculate initial period ratio/eccentricity just before tides are turned on.
    exit = 0
    inc = N+1+inner     #skip header to inner planet in resonance
    e0avg = np.zeros(0)
    p0avg = np.zeros(0)
    while exit != 1:
        temp2 = lines[inc]
        tempt = temp2.split("\t")
        tt = float(tempt[0])
        if tt > 70000 and tt < 80000:
            temp2out = lines[inc + (outer - inner)]
            temptout = temp2out.split("\t")
            e0avg = np.append(e0avg, float(tempt[2]))
            p0avg = np.append(p0avg, float(temptout[3])/float(tempt[3]))
        elif tt > 80000:
            e0 = np.append(e0, np.median(e0avg))
            p0 = np.append(e0, np.median(p0avg))
            exit = 1
        inc += N
    inp = np.mean(Pavg)
    simratio = outp/inp
    PRsim = np.append(PRsim,simratio)
    tempe = lines[length-N+inner]
    tempe = tempe.split("\t")
    time[k] = float(tempe[0])/1000000.
    e = float(tempe[2])
    #if e > 0.05:
    #    e = 0.05001
    e_in = np.append(e_in,e)   #eccentricity of inner planet (estimate of tides)

fos.close()
print 'the median eccentricity of Kepler planets pre-tides is ', np.median(e0)
#plot cdf
if early == 1:
    simp = p0
    sime = e0
    xmax = 0.16
else:
    simp = PRsim
    sime = e_in
    xmax = 0.06

binwidth = 0.0001
plt.hist(PRobs - 2.0, color='black', linewidth=2, bins=np.arange(0., 0.06 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label='Observed $\Delta$ from 2:1 MMR')
plt.hist(simp - 2.0, color='green', linewidth=2, bins=np.arange(0., 0.06 + binwidth, binwidth), histtype='step',cumulative='true', normed='true', label='Simulated $\Delta $ from 2:1 MMR after '+str(round(max(time)))+'Myr')
plt.hist(sime, color='orange', linewidth=2, bins=np.arange(0., 0.15 + binwidth, binwidth), histtype='step', cumulative='true', normed='true', label='Simulated e$_{inner}$')

plt.ylim([0,1.25])
plt.xlim([0,xmax])
plt.xlabel('P$_2$/P$_1$ - 2 or e')
plt.ylabel('counts, total='+str(len(PRsim)))
plt.title('Period ratios and Eccentricities of Planets Close to 2:1 MMR. '+ext)
plt.legend(loc='upper left',prop={'size':10})
plt.show()

#plot cdf the mathy way
#simp = np.sort(simp)
#sime = np.sort(sime)
#PRobs = np.sort(PRobs)
#yarr = np.arange(1,len(simp)+1, dtype=np.float)/len(simp)
#plt.plot(PRobs - 2.0, yarr, color='black', linewidth=2, label='Observed $\Delta$ from 2:1 MMR')
#plt.plot(simp - 2.0, yarr, color='green', linewidth=2, label='Simulated $\Delta $ from 2:1 MMR after '+str(round(max(time)))+'Myr')
#plt.plot(sime, yarr, color='orange', linewidth=2, label='Simulated e$_{inner}$')

#plot histogram
#binwidth = 0.01
#plt.hist(PRobs - 2.0, color='black', bins=np.arange(0., 0.06 + binwidth, binwidth), label='Observed $\Delta$ from 2:1 MMR')
#plt.hist(PRsim - 2.0, color='green', alpha=0.7, bins=np.arange(0., 0.06 + binwidth, binwidth), label='Simulated $\Delta$ from 2:1 MMR after '+str(round(max(time)))+'Myr')
#plt.hist(PRobs - PRsim, color='yellow', alpha=0.4, bins=np.arange(0., 0.06 + binwidth, binwidth), label='$\Delta_{Obs}$ - $\Delta_{Sim}$ (system-by-system basis)')
#plt.hist(e_in, alpha=0.3, color='yellow', bins=np.arange(0., 0.06 + binwidth, binwidth), label='Simulated e$_{inner}$')


