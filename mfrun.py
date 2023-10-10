# Preamble, loading necessary packages
import numpy as np
from configs import *
# import HEWG_inc_vard
import os
import platform

for frq in frequencies:
# frq = 73.0
    print("Working on frequency "+ str(frq))
    file = open('freq_used.py', 'w')
    file.write("frq= " + str(frq))
    file.close()
    if (platform.system() == 'Linux'):
        # os.system("python3 HEWG_inc_var_depth.py")
        os.system("python3 HEWG_tot_var_depth.py")
    elif (platform.system() == 'Windows'):
        os.system("python HEWG_inc_var_depth.py")
        os.system("python HEWG_tot_var_depth.py")
    else:
        print('You are using Mac, you need to fill this in') # Don't have a Mac to test on. Need to fill.