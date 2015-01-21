import multiprocessing as mp
import os

name=['TESTP3','TESTP5','TESTP10','TESTP10J']
os.system('make')

def execute(sysname):
    os.system('./nbody '+str(sysname))

#This
if __name__== '__main__':
    pool = mp.Pool(processes=3)
    args=[name[i] for i in xrange(0,3)]
    pool.map(execute, args)
    pool.close()
    pool.join()