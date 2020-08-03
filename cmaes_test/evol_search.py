"""
Created on July 02 10:52:00 2020
@author: carlos enrique gutierrez
DOYA UNIT - OIST    
"""


import numpy as np
import random
import os
import json
import pickle
import time
import sys
#sys.path.append('/work/DoyaU/pycma/')
import cma
#import subprocess
from settings_cmaes import *
#from mpi4py import MPI
import dotracking #gt_objfx

def main(brain_id):

    my_path_source = my_path + brain_id+'/moo_exvivo/light/'
    my_path_results = my_path + brain_id+'/moo_output/light_cmaes/'
    
    if True: #rank==0:
        # create a folder structure
        if not os.path.exists(my_path_results):
            os.makedirs(my_path_results)
        if not os.path.exists(my_path_results+'tracking/'):
            os.makedirs(my_path_results+'tracking/')
        if not os.path.exists(my_path_results+'conn_csv/'):
            os.makedirs(my_path_results+'conn_csv/')
        if not os.path.exists(my_path_results+'tracking2std/'):
            os.makedirs(my_path_results+'tracking2std/')
        if not os.path.exists(my_path_results+'evol/'):
            os.makedirs(my_path_results+'evol/')
        if not os.path.exists(my_path_results+'config/'):
            os.makedirs(my_path_results+'config/')

        evol_id = 0 #3#0
    
        #### Evolutiuonary algorithm
        population_size = 3#8
        es = cma.CMAEvolutionStrategy(3 * [0.5],0.01,{'popsize':population_size})
        #es = cma.CMAEvolutionStrategy([0.338618649636,0.532604109512,0.322303331729,0.398675048513,0.503575559129,0.616976662892],0.001,{'popsize':population_size})
        my_grid = ['angle','cutoff','minlength']  
        my_scale = [90., 0.2, 4.0] 


        while not es.stop():
            evol_id+=1
            
            with open(my_path_results+'evol/'+'obj_results_'+str(evol_id)+'.txt', 'a') as file:
                file.write('#brain_id, evol, obj, TPR global, FPR global, param1, param2, param3, TNR global, Spearman corr global'+'\n')
            
            population = es.ask()
            pop_aux = []
            for pop_i in population:
                pop_aux_j = []
                for pop_j in pop_i:
                    pop_aux_j.append(pop_j)
                pop_aux.append(pop_aux_j)
            
            data = {'pop': pop_aux, 'my_grid': my_grid, 'my_scale':my_scale,'my_path_source':my_path_source,
                       'my_path_results':my_path_results,'evol_id':evol_id,'brain_id':brain_id}
            with open(my_path_results+'config/'+'config_'+brain_id + '_' + str(evol_id)+'.txt', 'w') as outfile:
                json.dump(data, outfile, indent=4, sort_keys=True)
            
            script = open(my_path_results+'evol/'+'run_'+str(evol_id)+'.sh','w')
            c1 = 'mpiexec -np '+str(population_size)+' python dotracking.py '+my_path_results+'config/'+'config_'+brain_id + '_' + str(evol_id)+'.txt '#+\
            script.writelines(c1)
            script.close()
            command = 'source ' +my_path_results+'evol/'+'run_'+str(evol_id)+'.sh'
            print('submitting job: ',command)
            os.system(command)

            readyToGo = False
            while not readyToGo:
                with open(my_path_results+'evol/'+'obj_results_'+str(evol_id)+'.txt') as file:
                    lines = file.readlines()
                if len(lines)-1 >= (population_size):
                    readyToGo = True
                else:
                    print (len(lines)-1," lines -- Wait")
                    time.sleep(60)           

            my_index, my_obj = [],[]
            for j in lines[1:]:
                k = j.split(',')
                if k[0]!='#brain_id':
                    my_index.append(int(k[2]))
                    my_obj.append(float(k[3]))
        
            print(my_index,my_obj)
            sort_my_index = np.argsort(np.array(my_index))
            aux = np.array(my_obj)[sort_my_index]

            obj_values = [float(np_float) for np_float in aux]
            es.tell(population,obj_values)
            es.disp()
            r = es.result
            with open(my_path_results+'evol.txt','a') as file:
                file.write(str(r.fbest)+' '+' '.join(str(e) for e in r.xbest)+'\n')

        optimum = es.result_pretty()[:2]
        with open(my_path_results+'evol.txt','a') as file:
            file.write(str(optimum[1])+' '+' '.join(str(e) for e in optimum[0])+'\n')

        print('the end')

if __name__ == '__main__':
    brain_id = sys.argv[1] #read parameters
    main(brain_id)