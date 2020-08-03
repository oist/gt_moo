"""
Created on November 09 16:40:40 2018
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

# Training data set
#BRAINS = ['R01_0070_CM1180F','R01_0072_CM1176F','R01_0078_CM1347F','R01_0071_CM1178F',  
#          'R01_0029_CM696F']#,
BRAINS = ['R01_0030_CM690F','R01_0054_CM1060F','R01_0034_CM521F',
          'R01_0039_CM703F','R01_0033_CM694F']

#BRAINS = ['R01_0030_CM690F','R01_0070_CM1180F']#,'R01_0072_CM1176F']

# Test data set
#BRAINS = ['R01_0026_CM692F','R01_0040_CM710M','R01_0048_CM1011F',
#         'R01_0043_CM628F','R01_0053_CM1061F','R01_0046_CM1023M']

# parameters exploration ranges
grid={}

grid['angle'] = [25.,55.] #[10, 90]
grid['cutoff'] = [0.05,0.5] #[0.01, 1.0]
grid['minlength'] = [1.0,10.0] #[1.0, 18.0]
#grid['select'] = [100000, 100000000000] #not sure about this one.


my_path = '/Users/carlosgutierrez/paper_gt/database/' #'/work/DoyaU/carlos/dwi_pip/database/' #'/Users/carlosgutierrez/paper_gt/database/'#'/work/DoyaU/carlos/dwi_pip/database/' #dwi database path
mask = 'auto-brain_map_org_mask.nii' #whole-brain mask
my_mode = 'test' #'optimize' #'test' will run 5 times,  #'optimize' will run optimization

# how many runs of Global tracking per brain for the same parameters (results are averaged)
n_runs = 1 #1 
my_process = 'sequential' #'by_job' #sequential: will run each individual sequentially, by_job: will create a job and sned parallel individuals.
