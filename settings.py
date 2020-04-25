"""
Created on November 09 16:40:40 2018
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

# Training data set
BRAINS = ['R01_0070_CM1180F','R01_0072_CM1176F','R01_0078_CM1347F','R01_0071_CM1178F',
          'R01_0029_CM696F','R01_0030_CM690F','R01_0054_CM1060F','R01_0034_CM521F',
          'R01_0039_CM703F','R01_0033_CM694F']

# Test data set
#BRAINS = ['R01_0026_CM692F','R01_0040_CM710M','R01_0048_CM1011F',
#          'R01_0043_CM628F','R01_0053_CM1061F','R01_0046_CM1023M']


# Options for jobs building (change this as needed)
SLURM_OPTIONS_1 = ['#!/bin/bash \n\n',
                    '#SBATCH --time=7-00:00:00 \n',
                    '#SBATCH --mem 4gb \n',
                    '#SBATCH -c 1 \n',
                    '#SBATCH --partition=postproc1 \n', 
                    '#SBATCH --input=none \n']
SLURM_OPTIONS_1_1 = ['module load matlab/R2015a \n',
                    'module load python/3.5.0 \n',
                    'export PATH=\"/work/DoyaU/mrtrix3/bin:/work/DoyaU/deap:${PATH}\" \n',
                    'module load gcc/6.2.0 \n',
                    'module load eigen/3.3.1 \n',
                    'module load zlib/1.2.8 \n']
COMMAND_1 = 'sbatch ' 

SLURM_OPTIONS_2 = ['#!/bin/bash \n\n',
                    '#SBATCH --time=7-00:00:00 \n',
                    '#SBATCH --mem 30gb \n',
                    '#SBATCH -c 2 \n',
                    '#SBATCH --partition=compute \n',
                    '#SBATCH --input=none \n']
SLURM_OPTIONS_2_1 = ['module load matlab/R2015a \n',
                    'module load python/3.5.0 \n',
                    'export PATH=\"/work/DoyaU/mrtrix3/bin:/work/DoyaU/deap:${PATH}\" \n',
                    'module load gcc/6.2.0 \n',
                    'module load eigen/3.3.1 \n',
                    'module load zlib/1.2.8 \n']
COMMAND_2 = 'matlab -r \"addpath  /work/DoyaU/marcogt/;'

# parameters exploration ranges
grid={}
grid['width']= [0.01,0.15]  # -> [0.01,0.10]  
grid['length']= [0.24,0.65] # -> [0.32,0.65] 
grid['weight']= [0.01,0.22] # -> [0.01,0.13]
grid['chemPot2']=[0.05,0.6] # -> [0.01,0.22]
grid['connlike']=[0.5,6.0]  # -> [0.1,3.0]

my_path = '/work/DoyaU/carlos/dwi_pip/database/' #dwi database path
mask = 'auto-brain_map_org_mask.nii' #whole-brain mask at standard brain space
lod = 'suggestSparse'
my_mode = 'optimize' #'test' will run 5 times,  #'optimize' will run optimization

# how many runs of Global tracking per brain for the same parameters (results are averaged)
n_runs = 2 #1 


