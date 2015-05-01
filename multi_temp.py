import multiprocessing as mp
import os
import sys
import numpy as np

name = ['','a','b','c']
length=len(name)

os.system('make')

def execute(sysname):
    Qpfac = 100
    if sysname == '':
        Qpfac = 10000
    os.system('./nbody TESTP5m'+sysname+' '+str(Qpfac)+' 0.8')

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    #args=[name[i] for i in xrange(0,length)]    #xrange alreday goes to length-1
    #args=[name[i] for i in xrange(0,3)]
    args=[name[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
