"""
Created on June 24 13:56:40 2020
@author: carlos enrique gutierrez  carlosengutierrez@gmail.com
DOYA UNIT - OIST
Okinawa Institute of Science and TEchnology
"""

###### MAIN PROCESS, (Send a process for each brain) ########

from l_settings import *
import l_moo_loop
from mpi4py import MPI
comm = MPI.COMM_WORLD

#send n process, 1 process per brain.
print('sent process for: ',BRAINS[comm.rank] , '    mode: ',my_mode)
l_moo_loop.main(BRAINS[comm.rank] ,my_mode)






