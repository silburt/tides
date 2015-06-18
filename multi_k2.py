#Purpose of this multi is to explore how k2 affect the equilibrium eccentricity.

import multiprocessing as mp
import os
import sys
import numpy as np

k2=[1,2,4,8,0.5,0.25,0.125]
length = len(k2)

os.system('make')

def execute(k2val):
    os.system('./rebound TESTP3 '+str(k2val))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[k2[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()

