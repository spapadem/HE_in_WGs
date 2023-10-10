import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from configs import *


[X,Y] = np.meshgrid(xg,yg)
I = np.zeros((Nx,Ny),dtype='complex')

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

    # Dm = Dc

    omega = 2.*np.pi*frequencies[i]
    NPM = np.floor((2 * Dm * frequencies[i])/c0).astype(int)
    kn = omega/c0
    m = np.linspace(1,NPM,NPM)
    VV = np.sqrt(2/Dm) *np.sin(np.pi*np.outer(y_a,m)/Dm)
    lambdam = (m*np.pi / Dm)**2
    betam = np.sqrt(kn*kn - lambdam)
    Dbinv = np.diag(betam)
    A = np.zeros((NPM,NPM))
    for m in range(NPM):
        A[m,m] = (b/Dm) - (1/((m+1)*np.pi))* np.sin(((m+1)*np.pi*b)/Dm)*np.cos((2*(m+1)*np.pi*Dm/2)/Dm)
        for n in range(m+1,NPM):
            A[m,n] = (2/((n-m+1)*np.pi))*np.sin(((n-m+1)*np.pi*b)/(2*Dm))*np.cos(((n-m+1)*np.pi*Dm/2)/Dm) - (2/((n+m+1)*np.pi))*np.sin(((n+m+1)*np.pi*b)/(2*Dm))*np.cos(((n+m+1)*np.pi*Dm/2)/Dm)
    
    A = A + np.transpose(A) - np.diag(np.diag(A))

   


    Ua,Sa,Va = np.linalg.svd(A)
    Sa = Sa/Sa[0]
    Dvinv = np.diag(1/Sa)
    SJ = Dvinv @ Ua @ np.transpose(VV)
    pproj  = np.transpose(VV) @ usc @ VV
    # pproj = Ua @ Shat @ np.transpose(Ua)

    # J = 3
    # U,S,V = np.linalg.svd(pproj)
    # pproj = S[J] * np.reshape(U[:,J],(NPM,1)) @ np.transpose(np.reshape(U[:,J],(NPM,1)))
    for i in range(NPM):
        Gp = h*  Dbinv @ np.transpose(VV) @ np.reshape(G,(Nr,Nx*Ny))

    for m  in range(NPM):
        for n in range(NPM):
            I = I + np.conj(pproj[m,n])*Gp[m,:]*Gp[n,:]

I = I/np.max(abs(I))
print(xg)
print(yg)
plt.imshow(abs(np.reshape(I,(Nx,Ny))), cmap='jet', extent = [xg[0], xg[-1], yg[0], yg[-1]], aspect = 'equal',vmin = 0.4,vmax=1)
plt.colorbar()
plt.gca().add_patch(plt.Circle((x_sc,y_sc),b,fill=False,color='w'))
plt.savefig("image2.eps",format='eps')
plt.show()
