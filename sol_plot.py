import meshio
from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np


# msh = meshio.read("Mesh_test.msh")
# print(msh)
# x = msh.points[:,0]
# y = msh.points[:,1]


u = loadmat("data/data_incident")
u_inc = u['u']
u = loadmat("data/data_total")
u_tot = u['u']

W = 500
D = 200
PML_size = 20
x_min = 0 - PML_size
x_max = W + PML_size
y_min = 0
y_max = D
ax = plt.subplot2grid((3, 3), (0, 0))
ax.imshow(u_inc.real, cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Incident (real)')
ax = plt.subplot2grid((3, 3), (1, 0))
ax.imshow(u_inc.imag, cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Incident (imag)')
ax = plt.subplot2grid((3, 3), (2, 0))
ax.imshow(abs(u_inc), cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Incident (modulus)')
ax = plt.subplot2grid((3, 3), (0, 1))
ax.imshow(u_tot.real, cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Total (real)')
ax = plt.subplot2grid((3, 3), (1, 1))
ax.imshow(u_tot.imag, cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Total (imag))')
ax = plt.subplot2grid((3, 3), (2, 1))
ax.imshow(abs(u_tot), cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Total (modulus)')
ax = plt.subplot2grid((3, 3), (0, 2))
ax.imshow(u_tot.real-u_inc.real, cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Scattered (real)')
ax = plt.subplot2grid((3, 3), (1, 2))
ax.imshow(u_tot.imag-u_inc.imag, cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Scattered (imag))')
ax = plt.subplot2grid((3, 3), (2, 2))
ax.imshow(abs(u_tot-u_inc), cmap='jet', extent=[x_min, x_max, y_min, y_max], aspect='equal')
plt.title(f'Scattered (modulus)')
plt.show()


# ax = plt.subplot2grid((3, 1), (0, 0))
# ax.scatter(x,y, s= 1 , c=uf.real)  #extent=[x_discretization[0], x_discretization[-1], y_discretization[0], y_discretization[-1]], 
# ax = plt.subplot2grid((3, 1), (1, 0))
# ax.scatter(x,y, s= 1 , c=uf.imag)  #extent=[x_discretization[0], x_discretization[-1], y_discretization[0], y_discretization[-1]], 
# ax = plt.subplot2grid((3, 1), (2, 0))
# ax.scatter(x,y, s= 1 , c=np.abs(uf))  #extent=[x_discretization[0], x_discretization[-1], y_discretization[0], y_discretization[-1]], 
# plt.show()