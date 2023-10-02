# Preamble, loading necessary packages
import numpy as np
from configs import *
# import HEWG_inc_vard
import os

for frq in frequencies:
# frq = 73.0
    print("Working on frequency "+ str(frq))
    file = open('freq_used.py', 'w')
    file.write("frq= " + str(frq))
    file.close()
    os.system("python3 HEWG_inc_vard.py")
    os.system("python3 HEWG_tot_alt.py")

# for frq in frequencies:
#     print("Working on frequency "+ str(frq))
#     file = open('freq_used.py', 'w')
#     file.write("frq= "+str(frq))
#     file.close()
    # os.system("python HEWG_inc_vard.py")
 
