import multiprocessing as mp
import os
import sys
import numpy as np

name = ['a','b','c','d']
length=len(name)

os.system('make')

def execute(sysname):
    os.system('./nbody2 TESTP5m'+str(sysname)+' 100')

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    #args=[name[i] for i in xrange(0,length)]    #xrange alreday goes to length-1
    #args=[name[i] for i in xrange(0,3)]
    args=[name[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()