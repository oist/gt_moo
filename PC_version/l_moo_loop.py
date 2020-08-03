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
#sys.path.append('/apps/unit/DoyaU/deap/lib/python3.5/site-packages/')
from deap import base, tools, creator
from l_settings import *
import l_gt


def launch_job(params):
    if my_process == 'sequential':
        print('sending sequential tracking runs .... ')
        for index in np.arange(len(params['val'])):
            #param_val_i = params['val'][i]
            #param_individual_i = params['individual'][i]
            l_gt.main(params,n_runs,index)
        
    else:
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
    '''
    my_default = [45,0.1,2.0]
    init_l = list(np.random.normal([45,0.1,2.0],[0.01,0.001,0.01]))

    for i in np.arange(len(init_l)):
        if (init_l[i]<low[i] or init_l[i]>up[i]):
            init_l[i]=my_default[i]
    return init_l
    '''
    b = []
    lines = [line.split(' ') for line in open(my_path_results+'evol/moo_results.txt')]
    for k in np.arange(len(lines)):
        if lines[k][0]==str((gen1-1)):
            a = [float(lines[k][2][1:-1]),float(lines[k][3][:-1]),float(lines[k][4][:-2])] #angle, curoff, minlength
            if a[0]<25. or a[0]>55.: #angle
                a[0]=45.0
            if a[1]<0.05 or a[1]>0.5: #cutoff
                a[1] = 0.05
            if a[2]<1.0 or a[2]>10.:
                a[2]=2.0 
            b.append(a)
            print('aaaaaaa: ',a)
    return b[random.randint(0,7)]


def do_tracking(my_pop,group,brain_id):
    global my_count
    #my_count+=1
    my_count_list = []
    individual_list = []
    for i in my_pop:
        my_count+=1
        my_count_list.append(my_count)
        individual_list.append(list(i))
    
    # create a txt file with the parameters.
    params = {}
    params['individual']= individual_list #check this !!! #list(individual) 
    params['brain_id']=brain_id
    params['val'] = [brain_id+'_'+str(i)+'_'+str(gen1) for i in my_count_list]  #id list #brain_id+'_'+str(my_count)+'_'+str(gen1) #id
    params['my_path_source'] = my_path_source
    params['my_path_results'] = my_path_results
    params['mask']=mask
    params['group']=group
    params['gen']=str(gen1)

    with open(my_path_results+'config/'+'config_'+params['val'][0]+'_'+params['val'][-1]+'.txt', 'w') as outfile:
        json.dump(params, outfile, indent=4, sort_keys=True)
    
    # create and submit a job
    print('launching job for individual: ',params['val'])
    launch_job(params)
    return 'config_'+params['val'][0]+'_'+params['val'][-1]+'.txt'



def get_fitnesses(my_pop,group,brain_id):
    with open(my_path_results+'evol/'+'moo_results_'+str(gen1)+'_'+group+'.txt','a') as file: #create a summary for each fitness processing call
        file.write('#brain_id,TPR local,FPR local,TPR global,FPR global,param1,param2,param3,TNR local,Spearman corr local,TNR global,Spearman corr global' +'\n')

        
    a = []
    #send 1 job with many processes (1 process per individual)
    config_file = do_tracking(my_pop,group,brain_id)

    if my_process != 'sequential':
        readyToGo = False
        while not readyToGo: #wait until completion of fitness calculation for all the individuals
            with open(my_path_results+'evol/'+'moo_results_'+str(gen1)+'_'+group+'.txt','r') as file:
                lines = file.readlines()
            if len(lines)-1 == len(my_pop):
                readyToGo = True
            else:
                print("Fitness not ready, jobs completed: "+str(len(lines)-1)+" -- Please Wait")
                time.sleep(60)
    else:
        with open(my_path_results+'evol/'+'moo_results_'+str(gen1)+'_'+group+'.txt','r') as file:
            lines = file.readlines()

    for j in my_pop: #delivering fitnesses
        for k in np.arange(1,len(lines)):
            line_k = lines[k].split(',')
            if line_k[0]!='#brain_id':
                with open(my_path_results+'config/'+config_file) as json_file: #get the parameters values (an individual)
                    my_param = json.load(json_file)
                for q in np.arange(len(my_param['individual'])):
                    print('comparing  ',my_param['individual'][q], '  with  : ',j) #to retrieve the corresponding fitness
                    if my_param['individual'][q]==list(j):
                        print('equal!!!')
                        a.append(((1.0-float(line_k[1]),float(line_k[2]),1.0-float(line_k[3]),float(line_k[4])),line_k[0])) #f1,f2,f3,f4,id (objectives) (1-TPR_local,FPR_local,1-TRP_global,FPR_global,individual-id)
                        ##brain_id, TPR local, FPR local, TPR global, FPR global, param1, param2, param3, TNR local, Spearman corr local, TNR global, Spearman corr global,.
                        break
    return a

def optimize(grid, brain_id):
    global gen1

    ############## optimization settings ################
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,-1.0,-1.0,-1.0)) #weights=(-1.0,-1.0,-1.0,-1.0))  # change here depending on number of objectives. (minimization)
    creator.create("Individual", array.array, typecode='d', fitness=creator.FitnessMin)
    toolbox = base.Toolbox()
    BOUND_LOW = [grid['angle'][0],grid['cutoff'][0],grid['minlength'][0]]
    BOUND_UP = [grid['angle'][-1],grid['cutoff'][-1],grid['minlength'][-1]]

    NDIM = len(grid.keys()) # number of parameters.
    toolbox.register("evaluate", do_tracking)
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
    fitnesses = get_fitnesses(invalid_ind,'initial',brain_id) #process the fitness
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.__dict__ = {'id':fit[1]}
        ind.fitness.values = fit[0]
    pop = toolbox.select(pop, len(pop)) # This is just for assigning the crowding distance to the individuals (no selection is done)
    
    if gen1==0:
        with open(my_path_results+'evol/'+'moo_results.txt', 'a') as file:
            for ind in pop:
                file.write(str(gen1)+' '+str(ind)+' '+str(ind.fitness.values)+' '+str(ind.fitness.__dict__['id'])+' \n')

    for gen in np.arange(16,NGEN):#(9,NGEN):#(1,NGEN):
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
        fitnesses = get_fitnesses(invalid_ind,'first',brain_id) #process the fitness
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
                                with open(my_path+k+'/moo_output/light/evol/'+str(gen)+'.pkl','rb') as input:
                                    bests = pickle.load(input)
                            else:
                                with open(my_path+k+'/moo_output/light/evol/'+str(gen)+'.pkl','rb') as input:
                                    bests = bests + pickle.load(input) #accumale champions
                except:
                    bests=[]
                    pass
                print('gen: ',gen,'   gen1: ',gen1)
                print('bests: ',bests, 'and len: ',len(bests))
                if len(bests)==(MU_best*(len(BRAINS)-1)): #(champions == number of brains - 1)
                    readyToGo = True
                else:  
                    print(len(bests)," champions ready -- Wait until completion")
                    time.sleep(60) 
            print('MU_best: ',MU_best)
            #### all champions ready, the process continue ####   
            print('pop: ',pop)
            explorers = tools.selTournamentDCD(pop, MU_best*4) # Selection 3 #selTournamentDCD: number of individuals to select must be a multiple of 4
            print('explorers: ',explorers)
            explorers = explorers + bests # some explorers are matched to champions
            print('bests: ',bests)
            print('explorers + bests: ',explorers)
            for ind3, ind4 in zip(explorers[::2], explorers[1::2]):
                toolbox.mate(ind3, ind4) #matching
                del ind3.fitness.values, ind4.fitness.values
            
            explorers = explorers + bests # adding the original champions to the matched set
            invalid_ind = [ind for ind in explorers]
            fitnesses = get_fitnesses(invalid_ind,'second',brain_id) #process the fitness
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.__dict__ = {'id':fit[1]}
                ind.fitness.values = fit[0]
            pop = toolbox.select(pop+explorers,MU) #Selection 4 (from the whole set we select the best for the nex generation)

        with open(my_path_results+'evol/'+'moo_results.txt', 'a') as file: #record the evolution results
            for ind in pop:
                file.write(str(gen1)+' '+str(ind)+' '+str(ind.fitness.values)+' '+str(ind.fitness.__dict__['id'])+' \n')

def run_ft_with_fixed_params(test_params,brain_id,number_runs):
    global my_count   
    my_count_list = []
    individual_list = []
    for i in np.arange(number_runs):
        my_count+=1
        my_count_list.append(my_count)
        individual_list.append(test_params[i])
    
    # create a txt file with the parameters.
    params = {}
    params['individual']= individual_list #check this !!! #list(individual) 
    params['brain_id']=brain_id
    params['val'] = [brain_id+'_'+str(i)+'_'+str(gen1) for i in my_count_list]  #id list #brain_id+'_'+str(my_count)+'_'+str(gen1) #id
    params['my_path_source'] = my_path_source
    params['my_path_results'] = my_path_results
    params['mask']=mask
    params['group']='test'
    params['gen']=str(gen1)

    with open(my_path_results+'config/'+'config_'+params['val'][0]+'_'+params['val'][-1]+'.txt', 'w') as outfile:
        json.dump(params, outfile, indent=4, sort_keys=True)
    
    # create and submit a job
    print('launching job for individual: ',params['val'])
    launch_job(params)
    return 'config_'+params['val'][0]+'_'+params['val'][-1]+'.txt'



def main(brain_id,my_mode ):#if __name__ == '__main__':    
    #my_mode: def or opt  (default or optimized, used with fixed parameters.)
    global my_count# counter for job id, it must be re-initialized to a new value in case of a re-run of the optimization.
    my_count = 0
    global gen1 # evolution id.
    gen1 = 0 #16 #9 #8 #0
    global my_path_source # source data for each brain
    my_path_source = my_path + brain_id+'/moo_exvivo/light/'
    global my_path_results #'/moo_output_test_known_opt/'# # optimization results
    if my_mode=='optimize':
        my_path_results = my_path + brain_id+'/moo_output/light/'
    else:
        my_path_results = my_path + brain_id+'/moo_output/light_test_on_training_opt_1/'  #def or opt  ( for test training /light_test_on_test_def)
    
    # create folders to store the results:
    if not os.path.exists(my_path_results):
        os.makedirs(my_path_results)
    if not os.path.exists(my_path_results+'tracking/'): #full set of fibers at DWI space
        os.makedirs(my_path_results+'tracking/')
    if not os.path.exists(my_path_results+'tckinj/'): #sub-set of fibers at standard brain space (intersection of full set of fibers with injection point)
        os.makedirs(my_path_results+'tckinj/')
    if not os.path.exists(my_path_results+'evol/'): 
        os.makedirs(my_path_results+'evol/')
    if not os.path.exists(my_path_results+'tracking2std/'): #full set of fibers to standard brain space
        os.makedirs(my_path_results+'tracking2std/')
    if not os.path.exists(my_path_results+'trackingmapstd/'): #full set of fibers to standard brain space density maps
        os.makedirs(my_path_results+'trackingmapstd/')
    if not os.path.exists(my_path_results+'conn_csv/'): #connection matrices (full set of fibers at standard brain space)
        os.makedirs(my_path_results+'conn_csv/')
    if not os.path.exists(my_path_results+'config/'): #jobs parameters
        os.makedirs(my_path_results+'config/')
    if not os.path.exists(my_path_results+'jobs/'): #jobs files
        os.makedirs(my_path_results+'jobs/')
    if not os.path.exists(my_path_results+'map/'): #jobs files
        os.makedirs(my_path_results+'map/')

    with open(my_path_results+'evol/'+'moo_results.txt', 'a') as file:
        file.write('#gen ind f1 f2 f3 f4'+'\n')
    with open(my_path_results+'evol/'+'bests.txt', 'a') as file:
        file.write('#gen champions'+'\n')
    
    if my_mode=='optimize':
        optimize(grid,brain_id)
    if my_mode=='test': 
        #for i in range(5): #run 5 times for each brain.
        #angle:  38.15109891552099 5.945226453676253
        #cutoff:  0.03765836060482289 0.0176659492702961
        #minlength:  2.123194407929269 0.5688634266707114
        
        #angle:  32.18438616151345 6.286521520704324
        #cutoff:  0.046538331720147236 0.01185133384102836
        #minlength:  4.788743544615779 2.5184043117873243


        my_runs = 5
        test_params = []
        for i in np.arange(my_runs):
            #test_params.append([]) # for the case of default
            #test_params.append(list(np.random.normal([38.15109891552099,0.03765836060482289,2.123194407929269],[5.945226453676253,0.0176659492702961,0.5688634266707114]))) # optimized params.
            test_params.append(list(np.random.normal([32.18438616151345,0.046538331720147236,4.788743544615779],[6.286521520704324,0.01185133384102836,2.5184043117873243]))) # optimized params.

        run_ft_with_fixed_params(test_params,brain_id,my_runs)
print('the-end')
