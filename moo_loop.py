"""
Created on November 09 16:40:40 2018
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

import numpy as np
import random
import os
import json
import pickle
import time, array
import sys
sys.path.append('/apps/unit/DoyaU/deap/lib/python3.5/site-packages/')
from deap import base, tools, creator
from settings_okinawa import *


def launch_job(params):
    slurmOptions = SLURM_OPTIONS_2 + ['#SBATCH --job-name=job_'+params['val']+' \n',
                    '#SBATCH -e '+params['my_path_results']+'jobs/job_'+params['val']+'.out \n',
                    '#SBATCH -o '+params['my_path_results']+'jobs/job_'+params['val']+'.out \n'] +\
                    SLURM_OPTIONS_2_1 + ['time python3 gt.py '+params['my_path_results']+'config/config_'+params['val']+'.txt '+str(n_runs)+' \n']

    script = open(params['my_path_results']+'jobs/'+'job_'+params['val']+'.slurm','w')
    script.writelines(slurmOptions)
    script.close()
    command = COMMAND_1 + params['my_path_results']+'jobs/'+'job_'+params['val']+'.slurm'
    print('submitting job: ',command)
    os.system(command)


def uniform(low, up, size=None):
    my_default = [0.1,0.3,0.133,0.2,0.5]
    init_l = list(np.random.normal([0.1,0.3,0.133,0.2,0.5],[0.01,0.01,0.001,0.01,0.01])) # default.
    for i in np.arange(len(init_l)):
        if (init_l[i]<low[i] or init_l[i]>up[i]):
            init_l[i]=my_default[i]
    return init_l


def globaltracking(individual,group):
    global my_count
    my_count+=1

    # create a txt file with the parameters.
    params = {}
    params['individual']=list(individual)
    params['brain_id']=brain_id
    params['val'] = brain_id+'_'+str(my_count)+'_'+str(gen1) #id
    params['my_path_source'] = my_path_source
    params['my_path_results'] = my_path_results
    params['lod'] = lod
    params['mask']=mask
    params['group']=group
    params['gen']=gen1
    with open(my_path_results+'config/'+'config_'+params['val']+'.txt', 'w') as outfile:
        json.dump(params, outfile, indent=4, sort_keys=True)
    
    # create and submit a job
    if (os.path.isfile(params['my_path_results']+'jobs/job_'+params['val']+'.slurm')):
        print('job already done : ',params['val'])
    else:
        print('launching job for individual: ',params['val'])
        launch_job(params)


def get_fitnesses(my_pop,group):
    with open(my_path_results+'evol/'+'moo_results_'+str(gen1)+'_'+group+'.txt','a') as file: #create a summary for each fitness processing call
        file.write('#brain_id, TPR*, penalty, ratio, correlation, param1, param2, param3, param4, param5' +'\n')
        
    a = []
    for i in my_pop: #send q job for each individual
        globaltracking(i,group) 

    readyToGo = False
    while not readyToGo: #wait until completion of fitness calculation for all the individuals
        with open(my_path_results+'evol/'+'moo_results_'+str(gen1)+'_'+group+'.txt','r') as file:
            lines = file.readlines()
        if len(lines)-1 == len(my_pop):
            readyToGo = True
        else:
            print("Fitness not ready, jobs completed: "+str(len(lines)-1)+" -- Please Wait")
            time.sleep(60)

    for j in my_pop: #delivering fitnesses
        for k in np.arange(1,len(lines)):
            line_k = lines[k].split(',')
            if line_k[0]!='#brain_id':
                with open(my_path_results+'config/config_'+line_k[0]+'.txt') as json_file: #get the parameters values (an individual)
                    my_param = json.load(json_file)
                print('comparing  ',my_param['individual'], '  with  : ',j) #to retrieve the corresponding fitness
                if my_param['individual']==list(j):
                    a.append(((float(line_k[2]),-1.*float(line_k[3]),-1.*float(line_k[4]),-1.*float(line_k[1])),line_k[0])) #f1,f2,f3,f4,id (objectives) (penalty,ratio,right-side-correlation,weighted-TPR,individual-id)
                    break
    return a


def optimize(grid):
    global gen1

    ############## optimization settings ################
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,-1.0,-1.0,-1.0))  # change here depending on number of objectives. (minimization)
    creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)
    toolbox = base.Toolbox()
    BOUND_LOW = [grid['width'][0],grid['length'][0],grid['weight'][0],grid['chemPot2'][0],grid['connlike'][0]]
    BOUND_UP = [grid['width'][-1],grid['length'][-1],grid['weight'][-1],grid['chemPot2'][-1],grid['connlike'][-1]]
    NDIM = len(grid.keys()) # number of parameters.
    toolbox.register("evaluate", globaltracking)
    toolbox.register("attr_float", uniform, BOUND_LOW, BOUND_UP, NDIM)
    toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.attr_float)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxSimulatedBinaryBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0)
    toolbox.register("mutate", tools.mutPolynomialBounded, low=BOUND_LOW, up=BOUND_UP, eta=20.0, indpb=1.0/NDIM)
    toolbox.register("select", tools.selNSGA2)
    
    NGEN = 100  # number of generations
    MU = 8 #population size
    MU_best = 1 # number of champions to export to other brains optimizations
    CXPB = 0.2 # crossover probability
    pop = toolbox.population(n=MU) # generate the initial population
    
    invalid_ind = [ind for ind in pop if not ind.fitness.valid]
    fitnesses = get_fitnesses(invalid_ind,'initial') #process the fitness
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.__dict__ = {'id':fit[1]}
        ind.fitness.values = fit[0]
    pop = toolbox.select(pop, len(pop)) # This is just for assigning the crowding distance to the individuals (no selection is done)
    
    if gen1==0:
        with open(my_path_results+'evol/'+'moo_results.txt', 'a') as file:
            for ind in pop:
                file.write(str(gen1)+' '+str(ind)+' '+str(ind.fitness.values)+' '+str(ind.fitness.__dict__['id'])+' \n')

    for gen in range(NGEN):
        ######### generational process #############
        gen1=gen
        offspring = tools.selTournamentDCD(pop, len(pop))  # Tournament (Selection 1)
        offspring = [toolbox.clone(ind) for ind in offspring] 
        for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
            if random.random() <= CXPB:    #matching 
                toolbox.mate(ind1, ind2)
            toolbox.mutate(ind1)           #mutation
            toolbox.mutate(ind2)           #mutation
            del ind1.fitness.values, ind2.fitness.values
   
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = get_fitnesses(invalid_ind,'first') #process the fitness
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.__dict__ = {'id':fit[1]}
            ind.fitness.values = fit[0]
            
        pop = toolbox.select(pop + offspring, MU) # Selection 2
        pop_to_export = toolbox.select(pop,MU_best) # Selection of champion(s)
        with open(my_path_results+'evol/'+str(gen)+'.pkl','wb') as output:
            pickle.dump(pop_to_export,output,pickle.HIGHEST_PROTOCOL) #save champion(s) as pickle object
        with open(my_path_results+'evol/'+'bests.txt', 'a') as file: 
            file.write(str(gen)+' '+str(pop_to_export)+'\n')

        if gen > 1: # champions sharing starts from evolution 2 (it can be changed to start from begining of the generational process)
            readyToGo = False
            #retrieve champions from all the brains
            while not readyToGo:
                my_flag = True
                try:
                    for k in BRAINS:
                        if k!=brain_id:
                            if my_flag:
                                my_flag=False
                                with open(my_path+k+'/moo_output/evol/'+str(gen)+'.pkl','rb') as input:
                                    bests = pickle.load(input)
                            else:
                                with open(my_path+k+'/moo_output/evol/'+str(gen)+'.pkl','rb') as input:
                                    bests = bests + pickle.load(input) #accumale champions
                except:
                    bests=[]
                    pass
                if len(bests)==(MU_best*(len(BRAINS)-1)): #(champions == number of brains - 1)
                    readyToGo = True
                else:  
                    print(len(bests)," champions ready -- Wait until completion")
                    time.sleep(60) 
            
            #### all champions ready, the process continue ####   
            explorers = tools.selTournamentDCD(pop, MU_best*3) # Selection 3
            explorers = explorers + bests # some explorers are matched to champions
            for ind3, ind4 in zip(explorers[::2], explorers[1::2]):
                toolbox.mate(ind3, ind4) #matching
                del ind3.fitness.values, ind4.fitness.values
            
            explorers = explorers + bests # adding the original champions to the matched set
            invalid_ind = [ind for ind in explorers]
            fitnesses = get_fitnesses(invalid_ind,'second') #process the fitness
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.__dict__ = {'id':fit[1]}
                ind.fitness.values = fit[0]
            pop = toolbox.select(pop+explorers,MU) #Selection 4 (from the whole set we select the best for the nex generation)

        with open(my_path_results+'evol/'+'moo_results.txt', 'a') as file: #record the evolution results
            for ind in pop:
                file.write(str(gen1)+' '+str(ind)+' '+str(ind.fitness.values)+' '+str(ind.fitness.__dict__['id'])+' \n')

if __name__ == '__main__':    
    brain_id = sys.argv[1]
    my_count = 0 # counter for job id, it must be re-initialized to a new value in case of a re-run of the optimization.
    gen1 = 0 # evolution id.
    my_path_source = my_path + brain_id+'/moo_exvivo/' # source data for each brain
    my_path_results = my_path + brain_id+'/moo_output/' # optimization results
    
    # create folders to store the results:
    if not os.path.exists(my_path_results):
        os.makedirs(my_path_results)
    if not os.path.exists(my_path_results+'config/'): #jobs parameters
        os.makedirs(my_path_results+'config/')
    if not os.path.exists(my_path_results+'jobs/'): #jobs files
        os.makedirs(my_path_results+'jobs/')
    if not os.path.exists(my_path_results+'tracking/'): #full set of fibers at DWI space
        os.makedirs(my_path_results+'tracking/')
    if not os.path.exists(my_path_results+'tckinj/'): #sub-set of fibers at standard brain space (intersection of full set of fibers with injection point)
        os.makedirs(my_path_results+'tckinj/')
    if not os.path.exists(my_path_results+'map/'): #sub-set of fibers at standard brain space density maps
        os.makedirs(my_path_results+'map/')
    if not os.path.exists(my_path_results+'evol/'): #optimization results summary
        os.makedirs(my_path_results+'evol/')
    if not os.path.exists(my_path_results+'tracking2std/'): #full set of fibers to standard brain space
        os.makedirs(my_path_results+'tracking2std/')
    if not os.path.exists(my_path_results+'trackingmapstd/'): #full set of fibers to standard brain space density maps
        os.makedirs(my_path_results+'trackingmapstd/')
    if not os.path.exists(my_path_results+'conn_csv/'): #connection matrices (full set of fibers at standard brain space)
        os.makedirs(my_path_results+'conn_csv/')

    with open(my_path_results+'evol/'+'moo_results.txt', 'a') as file:
        file.write('#gen ind f1 f2 f3 f4'+'\n')
    with open(my_path_results+'evol/'+'bests.txt', 'a') as file:
        file.write('#gen champions'+'\n')
    
    optimize(grid)

print('the-end')
