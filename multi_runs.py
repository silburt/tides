#The purpose of this macro is to investigate the effect of different timesteps on the results of a system. It requires that the second argum ent in the "./nbody ..." sequence is the number of timesteps per orbit of innermost planet in problem.c, and make sure that it correctly links to readplanets.c. arg0 is now the name of the system

import multiprocessing as mp
import os
import sys
import numpy as np

runs=[('TESTP5m',0.5,20000),('TESTP5ma',0.2,100),('TESTP5mb',0.5,300),('TESTP5mc',0.5,500)]
length = len(runs)

os.system('make')

def execute(runvals):
    os.system('./rebound '+runvals[0]+' '+str(runvals[1])+' '+str(runvals[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
