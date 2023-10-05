import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from configs import *
import cmath


[X,Y] = np.meshgrid(x,y)
I = np.zeros(Nx*Ny,dtype='complex')
for  i in range(Nf):
    
    # Load incident field.
    ui = loadmat('data/inc_f'+ str(frequencies[i]) + '.mat')
    uinc = ui['u']
    # Load total field.
    ut = loadmat('data/tot_f'+ str(frequencies[i]) + '.mat')
    utot = ut['u']
   # Create scattered field.
    usc = utot - uinc

    Gf = loadmat('data/green_f' + str(frequencies[i]) + '.mat')
    G = Gf['u']

    omega = 2.*np.pi*frequencies[i]
    NPM = np.floor((2 * Dm * frequencies[i])/c0).astype(int)
    kn = omega/c0
    m = np.linspace(1,NPM,NPM)
    VV = np.sqrt(2/Dm) *np.sin(np.pi*np.outer(y_a,m)/Dm)
    lambdam = (m*np.pi / Dm)**2
    betam = np.sqrt(kn*kn - lambdam)
    Dbinv = np.diag(betam)
    B = h*(np.transpose(VV) @ VV)

    Ua,Sa,Va = np.linalg.svd(B)
    Dvinv = np.diag(1/Sa)
    SJ = Dvinv @ np.transpose(Ua) @ np.transpose(VV)
    Shat  = SJ @ usc @ np.transpose(SJ)
    pproj = Ua @ Shat @ np.transpose(Ua)
    Gp = h*np.transpose(VV) @ np.reshape(G,(Nr,Nx*Ny))

for m  in range(NPM):
    for n in range(NPM):
        I = I + np.conj(pproj[m,n])*Gp[m,:]*Gp[n,:]

plt.imshow(abs(np.reshape(I,(Nx,Ny))), cmap='jet', extent = [xg[0], xg[-1], yg[0], yg[-1]], aspect = 'equal')
plt.show()
print('Ok')