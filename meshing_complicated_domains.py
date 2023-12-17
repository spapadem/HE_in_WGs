# Preamble, loading necessary packages
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
from scipy.io import savemat
from scipy.interpolate import griddata
import numpy as np
import meshio
import ngsolve.internal as ngsint
# ngsint.viewoptions.drawoutline=0 # disable triangle outline when plotting.




Nbn = 151 # Number of nodes on the boundary.

# Parametric equation of the shape.
t = np.linspace(0,2*pi,Nbn)
t = t[0:-1]
n = 0 
rvals = np.abs(np.cos(n*t))**(np.sin(2*n*t))

# Calculate x and y (polar).
x = rvals*np.cos(t) 
y = rvals*np.sin(t)

# Creating the geometry.
geo = SplineGeometry()

# Create a list of points for the geometry.
P = np.zeros((t.shape[0],2))
P[:,0] = x
P[:,1] = y
pnts = [tuple(x) for x in P]
pp = [geo.AppendPoint(*pnt) for pnt in pnts]

# Creating edges.
curves = []
for i in range(Nbn-2):
    curves.append([["line",pp[i],pp[i+1]],"round"])
curves.append([["line",pp[-1],pp[0]],"round"])
[geo.Append(c,bc=bc) for c,bc in curves]

# Generate the mesh.
mesh = Mesh(geo.GenerateMesh(maxh=0.05))
mesh.Curve(3)
Draw(mesh)
print("Ok")
   
# Save the mesh as Gmsh2 format.
meshname = "Mesh_saving_test.msh"
mesh.ngmesh.Export(meshname,"Gmsh2 Format") # Saving the mesh file. Not needed if you choose to interpolate to a regular grid later on.