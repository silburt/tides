import sys
import matplotlib.pyplot as plt
import numpy as np
import math

rp = 0.0183*0.00464913       #Jupiter=0.1, Earth=0.01, Neptune = 0.035
mp = 0.00003   #Jupiter=.001, Earth=3e-6, Neptune = 5e-5
Ms = 1.         #Units of solar mass
a=0.057216999985679    #AU - THIS CHANGES OVER TIME!!!!
G = 1.          #AU^3/Msun*(yr/2pi)
Qp = 0.1/10.     #k_2/Q
pi = math.pi
com=Ms
R5a5 = rp*rp*rp*rp*rp/(a*a*a*a*a)
GM3a3 = (G*com*com*com/(a*a*a))**0.5
e=0.1
dt = 0.00782463923683634696

#time, a, e, i, Omega (long. of asc. node), omega, l (mean longitude), P, f
arg1=int(sys.argv[1])
arg2=int(sys.argv[2])
file_name=str(sys.argv[3])
names=['time (years)','de','dt','e','R5a5','GM3a3','a (AU)']
colors=['b','g','m','r','c']
fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=",")

N=1
for i in range(0,N): #range(0,N) only goes to N-1
    p=data[i::N]
    plt.plot(p[:,arg1], p[:,arg2], ''+colors[i])

lima=[0,0,0,0,0.0000000000000072692312315,73.0655024701,0.05721699999955]
limb=[0,0,0,0,0.0000000000000072692312325,73.0655024703,0.05721700000020]

if arg2 == 6:
    l = len(p[:,0])
    aa=np.zeros(l)
    for i in xrange(0,l):
        aa[i]=0.057217
    plt.plot(p[:,0],aa[:], 'k-.')

if arg2 ==1 or arg2==3:
    l = len(p[:,0])
    de_array=np.zeros(l)
    e_t=np.zeros(l)
    for i in xrange(0,l):
        a = p[i,6]
        R5a5 = rp*rp*rp*rp*rp/(a*a*a*a*a)
        GM3a3 = (G*com*com*com/(a*a*a))**0.5
        lne = -(9.*pi/2.)*Qp*GM3a3*R5a5/mp
        time = p[i,0]
        e_t[i] = math.e**(lne*time)*0.1 #the 0.1 is a constant. 
        de_array[i] = -dt*(9.*pi*0.5)*Qp*GM3a3*R5a5*e_t[i]/mp

if arg2==1:
    plt.plot(p[:,0],de_array[:], 'k-.')
if arg2==3:
    plt.plot(p[:,0],e_t[:], 'k-.')

if arg2 > 3:
    plt.ylim([lima[arg2],limb[arg2]])
plt.xlabel('' + names[arg1])
plt.ylabel('' + names[arg2])
plt.show()