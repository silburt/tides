#The purpose of this macro is to investigate the effect of different timesteps on the results of a system. It requires that the second argument in the "./nbody ..." sequence is the number of timesteps per orbit of innermost planet in problem.c, and make sure that it correctly links to readplanets.c. arg0 is now the name of the system

import multiprocessing as mp
import os
import sys
import numpy as np

arg0=str(sys.argv[1])
timesteps=[11,13,15,20,25]
length = len(timesteps)

os.system('make')

def execute(timefac):
    os.system('./nbody '+arg0+' '+str(timefac))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    #args=[name[i] for i in xrange(0,length)]    #xrange alreday goes to length-1
    #args=[name[i] for i in xrange(0,3)]
    args=[timesteps[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
