#The purpose of this macro is to output the Kepler systems in a nice table 1 format.

import sys
import numpy as np
import matplotlib.pyplot as plt
import math
import pylab

def blank(val):
    if val == 0:
        string = ''
    else:
        string = str(float("{0:.2f}".format(val)))
    return string

def blank2(val1):
    if val1 == '':
        string = ''
    else:
        string = '\pm'
    return string

ressystems=np.genfromtxt('../reso/2:1_systems.txt', dtype=("|S10")) #system names
ressystems = np.delete(ressystems,[1,7,8,16,21,29,32])
N_ressys = len(ressystems)

data=np.genfromtxt('../planets.txt', delimiter=',', dtype=(int,"|S10","|S10","|S10",int,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float), filling_values="0", skip_header=34)

output = open('Table1.txt','w')

RstoE = 109.17
nlines=data.shape[0]
count=0
j=0
while j < N_ressys:
    Np = data[count][4]
    if ressystems[j] == data[count][1]:
        index = count
        mp = data[index][19]
        mperr = data[index][20]
        Ms = data[index][15]
        Mserr = data[index][16]
        strmp = blank(mp)
        strmperr = blank(mperr)
        mppm = blank2(strmp)
        strMs = blank(Ms)
        if Ms > 0 and Mserr == 0:
            Mserr = 0.05
        strMserr = blank(Mserr)
        mspm = blank2(strMs)
        output.write(str(data[index][1])+' & '+str(data[index][2])+' & '+str(float("{0:.2f}".format(data[index][5])))+' & \checkmark & $'+strmp+' '+mppm+' '+strmperr+'$ & $'+str(float("{0:.2f}".format(data[index][22]*RstoE)))+' \pm '+str(float("{0:.2f}".format(data[index][23]*RstoE)))+'$ & '+str(data[index][4])+' & $'+strMs+' '+mspm+' '+strMserr+'$ & $'+str(data[index][17])+' \pm '+str(data[index][18])+'$ \\\  \n')
        for i in xrange(1,Np):
            index = i + count
            mp = data[index][19]
            mperr = data[index][20]
            strmp = blank(mp)
            strmperr = blank(mperr)
            mppm = blank2(strmp)
            output.write(' & '+str(data[index][2])+' & '+str(float("{0:.2f}".format(data[index][5])))+' & \checkmark & $'+strmp+' '+mppm+' '+strmperr+'$ & $'+str(float("{0:.2f}".format(data[index][22]*RstoE)))+' \pm '+str(float("{0:.2f}".format(data[index][23]*RstoE)))+'$ & & & \\\  \n')
        output.write('\hline \n')
        j += 1
    count += 1