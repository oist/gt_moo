#import numpy as np
#import random
#import os
import json
#import pickle
#import time
import sys
#sys.path.append('/work/DoyaU/pycma/')
#import cma
#import subprocess
#from settings_cmaes import *
from mpi4py import MPI
import gt_objfx




def main(data):
    
    comm = MPI.COMM_WORLD #send n process, 1 process per tracking.
    rank = comm.Get_rank()
    size = comm.Get_size()


    print('my rank is: ',rank)
    #data = {'pop': population[0], 'my_grid': my_grid, 'my_scale':my_scale,'my_path_source':my_path_source,
    #                   'my_path_results':my_path_results,'evol_id':evol_id,'rank':i,'brain_id':brain_id}
    params={}
    for i,k in enumerate(data['my_grid']):
        params[k] = data['pop'][rank][i] * data['my_scale'][i]

    params['my_path_source'] = data['my_path_source']
    params['my_path_results'] = data['my_path_results']
    params['evol_id']=str(data['evol_id'])
    params['val'] = data['brain_id'] + '_' + str(data['evol_id'])+'_'+str(rank)
    params['rank'] = str(rank)
    gt_objfx.main(params,1,1)

if __name__ == '__main__':
    data_file = sys.argv[1] #read parameters
    #data_path = sys.arg[2]

    with open(data_file) as json_file: #get the parameters values (an individual)
        data = json.load(json_file)

    main(data)


