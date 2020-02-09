"""
Created on November 09 16:40:40 2018
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

import os
from settings import *

###### MAIN PROCESS, (build and send a job for each brain) ########
for i in BRAINS:
    slurmOptions = SLURM_OPTIONS_1 + ['#SBATCH --job-name=job_'+i+' \n','#SBATCH -e '+i+'.out \n',
    '#SBATCH -o '+i+'.out \n'] + SLURM_OPTIONS_1_1 + ['time python3 moo_loop.py '+i+' '+my_mode+'\n']
    script = open(i+'.slurm','w')
    script.writelines(slurmOptions)
    script.close()
    os.system(COMMAND_1 +i+'.slurm')
    
print('seeds jobs submitted')
