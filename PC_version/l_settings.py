"""
Created on November 09 16:40:40 2018
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

# Training data set
BRAINS = ['R01_0070_CM1180F','R01_0072_CM1176F','R01_0078_CM1347F','R01_0071_CM1178F',  
          'R01_0029_CM696F']#,
#BRAINS = ['R01_0030_CM690F','R01_0054_CM1060F','R01_0034_CM521F',
#          'R01_0039_CM703F','R01_0033_CM694F']

# Test data set
#BRAINS = ['R01_0026_CM692F','R01_0040_CM710M','R01_0048_CM1011F','R01_0043_CM628F','R01_0053_CM1061F','R01_0046_CM1023M']

# parameters and exploration ranges
grid={}
grid['angle'] = [25.,55.] #[10, 90]
grid['cutoff'] = [0.05,0.5] #[0.01, 1.0]
grid['minlength'] = [1.0,10.0] #[1.0, 18.0]

my_mode = 'optimize' # will run the optimization. #'test' will run 5 times tractography with fixed parameters below (select default or optimized)
n_runs = 1 # how many runs of fiber-tracking per brain for the same parameters (results are averaged) (used in the optimization).
fixed_params = {'optimized':{'mean':[32.18438616151345,0.046538331720147236,4.788743544615779],
                              'std':[6.286521520704324,0.01185133384102836,2.5184043117873243]},
                'default':[],
                'mode':'optimized'}# or default (depending the case)
my_runs = 5 # used for the fixed parameters runs (no-optimization)

### Evolutionary algorithm parameters 
NGEN = 100  # max. number of generations
MU = 8 #population size
MU_best = 1 # number of champions to export to other brains optimizations
CXPB = 0.2 # crossover probability

##### Paths ######
my_path = '/Users/carlosgutierrez/paper_gt/database/'
source_path = '/moo_exvivo/light/'    # absolute path will be -->>>>  my_path + BRAINS[i] +  source_path
output_path =  '/moo_output/test/'    # absolute path will be -->>>>  my_path + BRAINS[i] +  output_path  (!!! CHANGE THIS TO AVOID OVERWRITTEN ON PREVIOUS RESULTS)



