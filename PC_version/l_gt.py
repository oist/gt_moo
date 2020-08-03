"""
Created on November 09 16:40:40 2018
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

import numpy as np
import subprocess
import sys
import SimpleITK as sitk
from l_settings import *
from scipy import stats
import scipy.io as sio
from numpy import matlib as mb
import json

def mapme(my_path_source,CM): #function for mapping to 20 x 104 tracer-based matrix (Code and tracer-based connectome from Skibbe H.)
    mapping = sio.loadmat(my_path_source+'atlas/'+'mat_mapping.mat')
    # the injection site regions we have in meso
    macro_srcs = np.squeeze(mapping['macro_srcs'])
    macro_srcs -= 1 # matlab2pyhton indexing
    #the mapping from high res atlas to low res atlas
    all_src_maps = np.squeeze(mapping['all_src_maps'])
    # number of unique targets - background
    valid_targets = np.unique(all_src_maps).shape[0]-1
    CM_macro = np.zeros((CM.shape[0],valid_targets))
    #run over all targets
    for b in range(1,valid_targets+1):    
        CM_macro[:,b-1] = np.sum(CM[:,all_src_maps == b],axis=1)    
    CMma = CM_macro[macro_srcs,:]
    CMma_norm = mb.repmat(np.sum(CMma,axis=1,keepdims=True),1,CMma.shape[1])
    CMma  = np.divide(CMma , CMma_norm)
    CMma_all = CM[macro_srcs,:]
    CMma_norm = mb.repmat(np.sum(CMma_all,axis=1,keepdims=True),1,CMma_all.shape[1])
    CMma_all  = np.divide(CMma_all , CMma_norm)
    return CMma, CMma_all


def maps_average(params,n_runs,index):
    my_maps,my_maps_1=[],[]
    for k in np.arange(int(n_runs)):
        my_temp = sitk.ReadImage(str(params['my_path_results']+'map/inj_map_'+params['val'][index]+'_'+str(k)+'.nii')) #average for fibers subset 
        my_temp_np = sitk.GetArrayFromImage(my_temp)
        my_maps.append(my_temp_np)

        my_temp_1 = sitk.ReadImage(str(params['my_path_results']+'trackingmapstd/map_'+params['val'][index]+'_'+str(k)+'.nii')) #average for fibers complete set
        my_temp_np_1 = sitk.GetArrayFromImage(my_temp_1)
        my_maps_1.append(my_temp_np_1)

    map_mean = np.mean(my_maps,axis=0)
    map_mean_nii = sitk.GetImageFromArray(map_mean)
    map_mean_nii.CopyInformation(my_temp)
    sitk.WriteImage(map_mean_nii,str(params['my_path_results']+'map/inj_map_'+params['val'][index]+'.nii'))
    del my_temp, my_temp_np, map_mean, map_mean_nii, my_maps
    map_mean_1 = np.mean(my_maps_1,axis=0)
    map_mean_nii_1 = sitk.GetImageFromArray(map_mean_1)
    map_mean_nii_1.CopyInformation(my_temp_1)
    sitk.WriteImage(map_mean_nii_1,str(params['my_path_results']+'trackingmapstd/map_'+params['val'][index]+'.nii'))
    del my_temp_1, my_temp_np_1, map_mean_1, map_mean_nii_1, my_maps_1


def main(params,n_runs,index):
    try: #if True: #try:
        for i in np.arange(int(n_runs)):
            # run tracking
            if len(params['individual'][index])==0:  # default case
                print('Runnning default case ......')
                c1 = 'tckgen '+  params['my_path_source'] + 'wmfod.mif ' + params['my_path_results'] + 'tracking/' + params['val'][index] +'_'+str(i)+'.tck' +\
                     ' -algorithm iFOD2 -seed_image ' + params['my_path_source']+'auto-brain_map_org_mask.nii' + ' -mask ' + params['my_path_source']+'auto-brain_map_org_mask.nii' +\
                     ' -grad ' + params['my_path_source']+'mrtrixbvecbval' +\
                     ' -select 300000 -force'
                print(c1)
                subprocess.call(c1, shell=True) # run tracking
            else:    
                c1 = 'tckgen '+  params['my_path_source'] + 'wmfod.mif ' + params['my_path_results'] + 'tracking/' + params['val'][index] +'_'+str(i)+'.tck' +\
                     ' -algorithm iFOD2 -seed_image ' + params['my_path_source']+'auto-brain_map_org_mask.nii' + ' -mask ' + params['my_path_source']+'auto-brain_map_org_mask.nii' +\
                     ' -grad ' + params['my_path_source']+'mrtrixbvecbval' +\
                     ' -cutoff ' + str(params['individual'][index][1]) + ' -minlength ' + str(params['individual'][index][2]) +\
                     ' -angle ' + str(params['individual'][index][0]) +\
                     ' -select 300000 -force'
                print(c1)
                subprocess.call(c1, shell=True)
            
            #move the complete set of fibers from dwi to standard brain space #tcktransform in latest mrtrix3 release
            #command = 'tcknormalise -force -info '+ params['my_path_results']+'tracking/'+params['val'][index]+'_'+str(i) + '.tck'+' '+\
            #my_path + params['brain_id']+ '/meta/tck_warp-dwi-to-std.mif ' + params['my_path_results'] + 'tracking2std/' + 'std_' + params['val'][index]+'_'+str(i) + '.tck'
            command = 'tcktransform -force -info '+ params['my_path_results']+'tracking/'+params['val'][index]+'_'+str(i) + '.tck'+' '+\
            my_path + params['brain_id']+ '/meta/tck_warp-dwi-to-std.mif ' + params['my_path_results'] + 'tracking2std/' + 'std_' + params['val'][index]+'_'+str(i) + '.tck'
            print('running ',command)
            subprocess.call(command, shell=True)
            
            #create connection matrix from the full set of fibers 
            TCK_FILE = params['my_path_results'] + 'tracking2std/' + 'std_' + params['val'][index]+'_'+str(i) + '.tck'
            MAP_NIFTI = my_path+'atlas/'+'mod-dir-Atlas_LabelMap_v2_1502.nii'
            CONNECTIVITY_CSV = params['my_path_results'] + 'conn_csv/' + 'std_' + params['val'][index]+'_'+str(i) + '.csv'
            command ="tck2connectome -symmetric -assignment_all_voxels " + TCK_FILE + " " + MAP_NIFTI + " " + CONNECTIVITY_CSV + ' -force'
            print('running ',command)
            subprocess.call(command, shell=True)

            #intersection with injection point in standard brain space
            command = 'tckedit -include '+params['my_path_source']+'inj_mask_std.nii '+params['my_path_results']+'tracking2std/'+'std_'+params['val'][index]+'_'+str(i)+'.tck'+' '+ \
            params['my_path_results']+'tckinj/inj_'+params['val'][index]+'_'+str(i)+'.tck -force'
            print('running: ',command)
            subprocess.call(command, shell=True)
            
            #subset of fibers to density map
            command = 'tckmap -template '+my_path+'model/mod-dir-average_brain/mod-dir-AverageBrain2_Isotropic.nii '+\
            params['my_path_results']+'tckinj/inj_'+ \
            params['val'][index]+'_'+str(i)+'.tck '+params['my_path_results']+'map/inj_map_'+params['val'][index]+'_'+str(i)+'.nii -force'
            print('running: ',command)
            subprocess.call(command, shell=True)
            
            #complete set of fibers to density map
            command = 'tckmap -template '+my_path+'model/mod-dir-average_brain/mod-dir-AverageBrain2_Isotropic.nii '+\
            params['my_path_results']+'tracking2std/'+ \
            'std_'+params['val'][index]+'_'+str(i)+'.tck '+params['my_path_results']+'trackingmapstd/map_'+params['val'][index]+'_'+str(i)+'.nii -force'
            print('running: ',command)
            subprocess.call(command, shell=True)

        # calculate and store average of density maps
        maps_average(params,n_runs,index)
        
        #load atlas
        img = sitk.ReadImage(my_path+'atlas/mod-dir-Atlas_LabelMap_v2_1502.nii')
        nda = sitk.GetArrayFromImage(img)
        vals = np.unique(nda)
        my_labels = np.loadtxt(my_path+'atlas/my_labels.txt')

        #load tracer
        my_boundary = 0.1 #remove noise from tracer
        my_tracer = sitk.ReadImage(str(params['my_path_source']+'tracer_masked_TC_org_2_MRI_std.nii.gz'))
        my_tracer_np = sitk.GetArrayFromImage(my_tracer)
        
        #load subset of fibers density map
        dfmri = sitk.ReadImage(str(params['my_path_results'] + 'map/inj_map_'+params['val'][index]+'.nii'))
        dfmri_np = sitk.GetArrayFromImage(dfmri)

        #create a single injection conectome (514 elements in 1 row) from the 1512 Rois Ids of the B/M atlas.
        connectome_tracer = []
        connectome_dfmri = []
        for i in my_labels: # construct rows based on atlas labels 
            idx = np.where(nda==i) #Roi coordinates
            element_value = len(my_tracer_np[idx][my_tracer_np[idx]>my_boundary]) #voxels counting (volume)
            #element_value = (my_tracer_np[idx][my_tracer_np[idx]>my_boundary]).sum() #voxels intensity sum
            connectome_tracer.append(element_value)
            element_value_dfmri = len(dfmri_np[idx][dfmri_np[idx]>0]) #voxels counting (volume)
            #element_value_dfmri = (dfmri_np[idx][dfmri_np[idx]>0]).sum() #voxels intensity sum
            connectome_dfmri.append(element_value_dfmri)
        del my_tracer
        del dfmri
        del my_tracer_np
        del dfmri_np

        # log transformation and normalization between 0 and 1
        connectome_dfmri_log = np.ma.log(connectome_dfmri) #log of the dti-matrix   
        connectome_dfmri_log = connectome_dfmri_log.filled(0.)
        connectome_dfmri_log_norm = (connectome_dfmri_log-np.min(connectome_dfmri_log))/(np.max(connectome_dfmri_log)-np.min(connectome_dfmri_log)) #normalization between 0-1 
        connectome_tracer_log = np.ma.log(connectome_tracer) #log of the tracer-matrix
        connectome_tracer_log = connectome_tracer_log.filled(0.)
        connectome_tracer_log_norm = (connectome_tracer_log-np.min(connectome_tracer_log))/(np.max(connectome_tracer_log)-np.min(connectome_tracer_log)) #normalization between 0-1 

        # Spearman correlation
        Spear_conn_row = stats.spearmanr(connectome_dfmri_log_norm.flatten(),connectome_tracer_log_norm.flatten()) #Spearman correlation coeficient

        #binarization of fibers and tracer: (preserving from 10% to 100% of the data)
        my_threshold = 0.1 
        connectome_tracer_log_norm_bin = np.array(list(connectome_tracer_log_norm))
        connectome_tracer_log_norm_bin[connectome_tracer_log_norm_bin<=my_threshold]=0.
        connectome_tracer_log_norm_bin[connectome_tracer_log_norm_bin>0]=1.0
        connectome_dfmri_log_norm_bin = np.array(list(connectome_dfmri_log_norm))
        connectome_dfmri_log_norm_bin[connectome_dfmri_log_norm_bin<=my_threshold]=0.
        connectome_dfmri_log_norm_bin[connectome_dfmri_log_norm_bin>0]=1.0

        # Individuaal TPR and FPR
        TPR_row = (np.logical_and(connectome_dfmri_log_norm_bin==1, connectome_tracer_log_norm_bin==1)).sum()/connectome_tracer_log_norm_bin.sum()
        TNR_row = (np.logical_and(connectome_dfmri_log_norm_bin==0, connectome_tracer_log_norm_bin==0)).sum()/len(connectome_tracer_log_norm_bin[connectome_tracer_log_norm_bin==0])
        FPR_row = 1. - TNR_row

        # Global TPR and FPR (as calculated in the previous version, Fig. 9 of the paper.)
        # load averaged DTI-based connection matrix (for the full set of fibers)
        aux =[]
        for k in np.arange(int(n_runs)):
            my_data = np.genfromtxt(params['my_path_results'] + 'conn_csv/' + 'std_' + params['val'][index]+'_'+str(k)+'.csv', delimiter=' ')
            aux.append(my_data)
        my_data_mean = np.mean(aux,axis=0)
        del my_data, aux

        #load tracer injections matrix.
        CMme_orig = sio.loadmat(my_path + 'atlas/'+'full_mat.mat')['CMme_full']#('/Users/carlosengutierrez/datathon/check_models/dwi_pip/atlas_Henrik/'+'CM_meso_norm.mat')['CMme']
        idx = np.where(np.sum(CMme_orig,axis=0)==0.)[0]
        CMme = np.delete(CMme_orig,idx,axis=1) 
        CMme = CMme**0.25   #transformation

        # Map dfrmi connectome to tracer matrix structure.
        CMma, CMma_all = mapme(my_path,my_data_mean)
        CMma = np.delete(CMma_all,idx,axis=1)
        CMma = np.delete(CMma,14,0)
        CMma = CMma**0.25  #transformation
        del my_data_mean
        
        #global TPR and FPR             
        my_threshold = 0.002 #0. #0.002 #0.1 #for removing some noise in matrices building.
        my_and = np.logical_and(CMme>my_threshold,CMma>my_threshold)
        my_and_2 = np.logical_and(CMme<=my_threshold,CMma<=my_threshold)

        TPR_global = my_and[my_and==1].shape[0]/float(CMme[CMme>my_threshold].shape[0])
        TNR_global = my_and_2[my_and_2==1].shape[0]/float(CMme[CMme<=my_threshold].shape[0])
        FPR_global = 1.-(my_and_2[my_and_2==1].shape[0]/float(CMme[CMme<=my_threshold].shape[0]))

        # Global Spearman correlation
        CMma_log = np.ma.log(CMma) #log of the tracks and normalize between 0-1         
        #CMma_log = CMma_log.filled(0.)  ##removing this, it seems it doesn't work for this case
        CMma_log_norm = (CMma_log-np.min(CMma_log))/(np.max(CMma_log)-np.min(CMma_log))
        CMme_log = np.ma.log(CMme) #log of the tracks.
        #CMme_log = CMme_log.filled(0.)  ##removing this, it seems it doesn't work for this case
        CMme_log_norm = (CMme_log-np.min(CMme_log))/(np.max(CMme_log)-np.min(CMme_log))
        Spear_conn_global = stats.spearmanr(CMma_log_norm.flatten(),CMme_log_norm.flatten())
        
        ############ Record the results ###############.
        if len(params['individual'][index])==0:
            params['individual'][index].append('test')
            params['individual'][index].append('test')
            params['individual'][index].append('test')
        

        with open(params['my_path_results']+'evol/'+'moo_results_'+str(params['gen'])+'_'+params['group']+'.txt','a') as file:
            #brain_id, TPR local, FPR local, TPR global, FPR global, param1, param2, param3, TNR local, Spearman corr local, TNR global, Spearman corr global,.....
            file.write(params['val'][index]+','+str(TPR_row)+','+str(FPR_row)+','+\
            str(TPR_global)+','+str(FPR_global)+','+\
            str(params['individual'][index][0])+','+str(params['individual'][index][1])+','+ str(params['individual'][index][2])+','+\
            str(TNR_row)+','+str(Spear_conn_row[0])+','+\
            str(TNR_global)+','+str(Spear_conn_global[0])+'\n') #additional data



        #### remove files ####
        #complete set of fibers to density map
        command = 'rm -r '+ params['my_path_results']+'tracking2std/'+ \
        'std_'+params['val'][index]+'_'+str(0)+'.tck'
        print('running: ',command)
        subprocess.call(command, shell=True)

        command = 'rm -r '+params['my_path_results']+'trackingmapstd/map_'+ \
        params['val'][index]+'_'+str(0)+'.nii'
        print('running: ',command)
        subprocess.call(command, shell=True)

    except:
    #else:
        print('tracking failed, going to the except branch')
        with open(params['my_path_results']+'evol/'+'moo_results_'+str(params['gen'])+'_'+params['group']+'.txt','a') as file:
            #brain_id, TPR local, FPR local, TPR global, FPR global, param1, param2, param3, TNR local, Spearman corr local, TNR global, Spearman corr global,.....
            file.write(params['val'][index]+','+str(0)+','+str(0)+','+\
            str(0)+','+str(0)+','+\
            str(params['individual'][index][0])+','+str(params['individual'][index][1])+','+ str(params['individual'][index][2])+','+\
            str(0)+','+str(0)+','+\
            str(0)+','+str(0)+'\n') #additional data
    
    

print('end')


        
