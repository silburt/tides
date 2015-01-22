import multiprocessing as mp
import os
import sys
import numpy as np

arg0=str(sys.argv[1])
arg1=str(sys.argv[2])

name=np.genfromtxt('reso/'+str(arg0)+':'+str(arg1)+'_systems.txt', dtype=("|S10"), delimiter=',')
length=name.shape[0]

os.system('make')

def execute(sysname):
    os.system('./nbody '+str(sysname)+' --tmax=100000 --tide_delay=80000')

#This is the main multiprocess execution
if __name__== '__main__':
    pool = mp.Pool(processes=3)
    #args=[name[i] for i in xrange(0,length)]    #xrange alreday goes to length-1
    args=[name[i] for i in xrange(0,3)]
    pool.map(execute, args)
    pool.close()
    pool.join()