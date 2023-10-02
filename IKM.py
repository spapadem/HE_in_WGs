import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from configs import *
t = loadmat('data/mesh_points_inc.mat')
p = t['u']

I = 0
i=0
for  i in range(1):
    # i
    # Load incident field.
    ui = loadmat('data/inc_f'+str(frequencies[i])+'.mat')
    uinc = ui['u']

    # Load total field.
    ut = loadmat('data/tot_f'+str(frequencies[i])+'.mat')
    utot = ut['u']

    # Create scattered field.
    usc = utot - uinc

    Gf = loadmat('data/green_f'+str(frequencies[i])+'.mat')
    G = Gf['u']

    for m  in range(usc.shape[0]):
        for n in range(usc.shape[1]):
            I = I + np.conjugate(usc[m,n])*(G[m,:])*(G[n,:])


plt.scatter(p[:,0],p[:,1],25,np.abs(I))
plt.show()

