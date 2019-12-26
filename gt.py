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
from settings import *
from scipy import stats
import scipy.io as sio
from numpy import matlib as mb

def mapme(my_path_source,CM): #function for mapping to 20 x 104 tracer-based matrix (Code from Skibbe H.)
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


def maps_average(params,n_runs):
    my_maps,my_maps_1=[],[]
    for k in np.arange(int(n_runs)):
        my_temp = sitk.ReadImage(str(params['my_path_results']+'map/inj_map_'+params['val']+'_'+str(k)+'.nii')) #average for fibers subset 
        my_temp_np = sitk.GetArrayFromImage(my_temp)
        my_maps.append(my_temp_np)

        my_temp_1 = sitk.ReadImage(str(params['my_path_results']+'trackingmapstd/map_'+params['val']+'_'+str(k)+'.nii')) #average for fibers complete set
        my_temp_np_1 = sitk.GetArrayFromImage(my_temp_1)
        my_maps_1.append(my_temp_np_1)

    map_mean = np.mean(my_maps,axis=0)
    map_mean_nii = sitk.GetImageFromArray(map_mean)
    map_mean_nii.CopyInformation(my_temp)
    sitk.WriteImage(map_mean_nii,str(params['my_path_results']+'map/inj_map_'+params['val']+'.nii'))
    del my_temp, my_temp_np, map_mean, map_mean_nii, my_maps
    map_mean_1 = np.mean(my_maps_1,axis=0)
    map_mean_nii_1 = sitk.GetImageFromArray(map_mean_1)
    map_mean_nii_1.CopyInformation(my_temp_1)
    sitk.WriteImage(map_mean_nii_1,str(params['my_path_results']+'trackingmapstd/map_'+params['val']+'.nii'))
    del my_temp_1, my_temp_np_1, map_mean_1, map_mean_nii_1, my_maps_1


def main(params,n_runs):
    try:
        for i in np.arange(int(n_runs)):
            # run global tracking
            c1 = COMMAND_2 + 'trackme(\'folder\',\''+params['my_path_source']+'\',\'lod\',\''+params['lod']+'\',\'file_mask\',\''+params['my_path_source']+params['mask']+'\',\'params\',{\'width\','+str(params['individual'][0])+ \
            ',\'length\','+str(params['individual'][1])+',\'weight\','+str(params['individual'][2])+ \
            ',\'chemPot2\','+ str(params['individual'][3])+',\'chemPot1\','+str(0.0)+',\'connlike\','+str(params['individual'][4])+ \
            '},\'ofile\',\''+params['my_path_results']+'tracking/'+params['val']+'_'+str(i)+'\');exit;\"'
            print(c1)
            subprocess.call(c1, shell=True)
            
            #move the complete set of fibers from dwi to standard brain space
            command = 'tcknormalise -force -info '+ params['my_path_results']+'tracking/'+params['val']+'_'+str(i) + '.tck'+' '+\
            params['my_path_source'][:-11]+ 'meta/tck_warp-dwi-to-std.mif ' + params['my_path_results'] + 'tracking2std/' + 'std_' + params['val']+'_'+str(i) + '.tck'
            print('running ',command)
            subprocess.call(command, shell=True)
            
            #create connection matrix from the full set of fibers 
            TCK_FILE = params['my_path_results'] + 'tracking2std/' + 'std_' + params['val']+'_'+str(i) + '.tck'
            MAP_NIFTI = my_path+'atlas/'+'mod-dir-Atlas_LabelMap_v2_1502.nii'
            CONNECTIVITY_CSV = params['my_path_results'] + 'conn_csv/' + 'std_' + params['val']+'_'+str(i) + '.csv'
            command ="tck2connectome -symmetric -assignment_all_voxels " + TCK_FILE + " " + MAP_NIFTI + " " + CONNECTIVITY_CSV
            print('running ',command)
            subprocess.call(command, shell=True)

            #intersection with injection point in standard brain space
            command = 'tckedit -include '+params['my_path_source']+'inj_mask_std.nii '+params['my_path_results']+'tracking2std/'+'std_'+params['val']+'_'+str(i)+'.tck'+' '+ \
            params['my_path_results']+'tckinj/inj_'+params['val']+'_'+str(i)+'.tck -force'
            print('running: ',command)
            subprocess.call(command, shell=True)
            
            #subset of fibers to density map
            command = 'tckmap -template '+my_path+'model/mod-dir-average_brain/mod-dir-AverageBrain2_Isotropic.nii '+\
            params['my_path_results']+'tckinj/inj_'+ \
            params['val']+'_'+str(i)+'.tck '+params['my_path_results']+'map/inj_map_'+params['val']+'_'+str(i)+'.nii -force'
            print('running: ',command)
            subprocess.call(command, shell=True)
            
            #complete set of fibers to density map
            command = 'tckmap -template '+my_path+'model/mod-dir-average_brain/mod-dir-AverageBrain2_Isotropic.nii '+\
            params['my_path_results']+'tracking2std/'+ \
            'std_'+params['val']+'_'+str(i)+'.tck '+params['my_path_results']+'trackingmapstd/map_'+params['val']+'_'+str(i)+'.nii -force'
            print('running: ',command)
            subprocess.call(command, shell=True)
        
        
        # calculate and store average of density maps
        maps_average(params,n_runs)
        
        ###### load data sources ############
        # load tracer
        my_boundary=0.1 #threshold to remove background
        template = sitk.ReadImage(str(params['my_path_source']+'tracer_masked_TC_org_2_MRI_std.nii.gz'))
        tracer_np = sitk.GetArrayFromImage(template)
        tracer_np[tracer_np>my_boundary]=99
        tracer_np[tracer_np<=my_boundary]=0
        tracer_np[tracer_np==99]=1 #image binarization (tracer image mask)
        tracer_np_weigthed = sitk.GetArrayFromImage(template) # original tracer image (weighted pixels)
        tracer_np_weigthed[tracer_np_weigthed<=my_boundary]=0 
        del template

        # load whole-brain mask
        img_mask = sitk.ReadImage(str(my_path+'mask/mask_std.nii')) 
        img_mask_np = sitk.GetArrayFromImage(img_mask)
        del img_mask
        
        #load subset of fibers density map
        tracks = sitk.ReadImage(str(params['my_path_results']+'map/inj_map_'+params['val']+'.nii'))
        tracks_np = sitk.GetArrayFromImage(tracks)
        tracks_np[tracks_np>0]=99
        tracks_np[tracks_np<=0]=0
        tracks_np[tracks_np==99]=1 # binarization
        del tracks


        ############# weighted TPR based on distances to the center of the inj. point and pixels intensity ################
        # center of mass of the injection point center
        img_inj = sitk.ReadImage(str(params['my_path_source']+'inj_center_TC_org_2_MRI_std.nii.gz'))
        img_inj_np = sitk.GetArrayFromImage(img_inj)
        img_inj_np[img_inj_np>0.1]=99
        img_inj_np[img_inj_np<=0.1]=0
        img_inj_np[img_inj_np==99]=1 # image binarization (injection center mask)
        roi_sum = np.sum(img_inj_np[img_inj_np==1])
        xyz = np.where(img_inj_np==1)
        cx,cy,cz = int(np.sum(xyz[0])/roi_sum), int(np.sum(xyz[1])/roi_sum), int(np.sum(xyz[2])/roi_sum) #center of mass coordinates.
        del img_inj, img_inj_np, roi_sum, xyz

        # distances from tracer voxels to center of mass
        xyz_tracer = np.where(np.logical_and(tracer_np==1,img_mask_np==1)==1)
        dist=np.sqrt((xyz_tracer[0]-cx)**2+(xyz_tracer[1]-cy)**2+(xyz_tracer[2]-cz)**2)
        dist_norm = dist/dist.max() #distances normalization

        # pixels weight from raw tracer intensity
        tracer_np_weigthed = tracer_np_weigthed/tracer_np_weigthed.max() #normalization
        
        # weights (pixels intensity and distance) combination for all NP (positives)
        dist_int_tracer = np.multiply(dist_norm,tracer_np_weigthed[xyz_tracer])
        dist_int_tracer = dist_int_tracer/np.sum(dist_int_tracer) #normalization

        # weighted TPR
        den_tra_eq = np.logical_and(tracks_np==tracer_np,tracer_np==1)
        den_tra_eq_neg = np.logical_and(tracks_np==tracer_np,tracer_np<1)
        xyz_den = np.where(np.logical_and(den_tra_eq,img_mask_np==1)==1)
        xyz_den_stack = np.vstack((np.vstack((xyz_den[0],xyz_den[1])),xyz_den[2])) # TP  (True positives)
        xyz_tracer_stack = np.vstack((np.vstack((xyz_tracer[0],xyz_tracer[1])),xyz_tracer[2])) # NP (all positives)
        xyz_den_stack = xyz_den_stack.T
        xyz_tracer_stack = xyz_tracer_stack.T
        del xyz_tracer, xyz_den
        cost_1 = 0
        for k in np.arange(len(xyz_den_stack)):
            cost_1+= dist_int_tracer[np.where(np.all(xyz_tracer_stack==xyz_den_stack[k,:],axis=1))[0][0]]  # sum of TP weights
        del xyz_den_stack, xyz_tracer_stack


        ############# Ratio ################
        # non-weighted TPR and FPR
        TPR = np.sum(np.logical_and(den_tra_eq,img_mask_np==1))/float(np.sum(np.logical_and(tracer_np==1,img_mask_np==1))) # True positive
        SPC = np.sum(np.logical_and(den_tra_eq_neg,img_mask_np==1))/float(np.sum(np.logical_and(tracer_np<=0,img_mask_np==1))) # True negative
        FPR = 1. - SPC # False positive
        del den_tra_eq_neg, den_tra_eq, img_mask_np, tracks_np, tracer_np

        mu_n = np.loadtxt('./my_epsilon.txt') #mean number of TN (training data set)
        mu_p = np.loadtxt('./my_gamma.txt')  #mean number of TP (training data set)
        tolerance = (0.01*mu_p*0.6)/mu_n #tolerance model
        cost_2 = cost_1/(FPR + tolerance) #ratio


        ############# correlation between DTI and tracer based connection matrices (right-side) ################
        # averaged DTI-based connection matrix (full set of fibers)
        aux =[]
        for k in np.arange(int(n_runs)):
            my_data = np.genfromtxt(params['my_path_results'] + 'conn_csv/' + 'std_' + params['val']+'_'+str(k)+'.csv', delimiter=' ')
            aux.append(my_data)
        my_data_mean = np.mean(aux,axis=0)
        del my_data, aux

        CMma, CMma_all = mapme(my_path,my_data_mean) #map connection matrix to 20x104 tracer connection matrix.
        CMme = sio.loadmat(my_path+'atlas/'+'CM_meso_norm.mat')['CMme'] # load tracer matrix
        del my_data_mean
        CMma_log = np.ma.log(CMma) #log of the dti-matrix   
        CMma_log_norm = (CMma_log-np.min(CMma_log))/(np.max(CMma_log)-np.min(CMma_log)) #normalization between 0-1 
        
        CMme_log = np.ma.log(CMme) #log of the tracer-matrix
        CMme_log_norm = (CMme_log-np.min(CMme_log))/(np.max(CMme_log)-np.min(CMme_log)) #normalization between 0-1 
        # Spearman correlation of the right hemisphere (projections)
        Spear_conn = stats.spearmanr(CMma_log_norm.flatten(),CMme_log_norm.flatten())
        cost_3 = stats.spearmanr(CMma_log_norm[:,:52].flatten(),CMme_log_norm[:,:52].flatten())

        
        ############# Penalty ################
        mask_penalty = sitk.ReadImage(str(my_path+'model/mask_penalty.nii')) #transversal mask
        mask_penalty_np = sitk.GetArrayFromImage(mask_penalty)
        del mask_penalty

        #load complete tracking density map 
        den_map = sitk.ReadImage(str(params['my_path_results']+'trackingmapstd/map_'+params['val']+'.nii'))
        den_map_np = sitk.GetArrayFromImage(den_map)
        del den_map

        #sum and normalization
        cost_4 = np.sum(np.logical_and(den_map_np>0,mask_penalty_np>0))/float(len(mask_penalty_np[mask_penalty_np>0])) #penalty
        del mask_penalty_np, den_map_np


        ############ Record the results ###############.
        with open(params['my_path_results']+'evol/'+'moo_results_'+str(params['gen'])+'_'+params['group']+'.txt','a') as file:
            #brain_id, TPR*, penalty, ratio, correlation, param1, param2, param3, param4, param5 (objectives and parameters)
            file.write(params['val']+','+str(cost_1)+','+str(cost_4)+','+str(cost_2)+','+str(cost_3[0])+','+\
            str(params['individual'][0])+','+str(params['individual'][1])+','+ str(params['individual'][2])+','+\
            str(params['individual'][3])+','+str(params['individual'][4])+'\n')            
    except:
        print('tracking failed, going to the except branch')
        with open(params['my_path_results']+'evol/'+'moo_results_'+str(params['gen'])+'_'+params['group']+'.txt','a') as file:
            #brain_id, TPR*, penalty, ratio, correlation, param1, param2, param3, param4, param5 (objectives and parameters)
            file.write(params['val']+','+str(0)+','+str(1)+','+str(0)+','+str(0)+','+\
            str(params['individual'][0])+','+str(params['individual'][1])+','+str(params['individual'][2])+','+\
            str(params['individual'][3])+','+str(params['individual'][4])+'\n')
        
print('owari')


if __name__ == '__main__':
    filename = sys.argv[1] #read parameters
    n_runs = sys.argv[2] #number of runs.
    with open(filename) as json_file:  
        params = json.load(json_file)
    main(params,n_runs)

        
