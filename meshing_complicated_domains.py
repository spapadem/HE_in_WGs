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


# Creating the waveguide geometry.
geo = SplineGeometry()

Nbn = 151
n = 1.5
t = np.linspace(0,2*pi,Nbn)
t = t[0:-1]
rvals = np.abs(np.cos(n*t))**(np.sin(2*n*t))

x = rvals*np.cos(t)
y = rvals*np.sin(t)

P = np.zeros((t.shape[0],2))
P[:,0] = x
P[:,1] = y
pnts = [tuple(x) for x in P]

pp = [geo.AppendPoint(*pnt) for pnt in pnts]
curves = []
for i in range(Nbn-2):
    curves.append([["line",pp[i],pp[i+1]],"round"])
curves.append([["line",pp[-1],pp[0]],"round"])
[geo.Append(c,bc=bc) for c,bc in curves]
# pnts = [pnt for pnt in (x,y)]

# geo2 = SplineGeometry()
# p1,p2,p3,p4,p5,p6,p7,p8,p9 = [geo2.AppendPoint(*pnt) for pnt in pnts]
# curves2 = [[["line",p1,p2],"top"],
#           [["line",p2,p3],"right"],
#           [["line",p3,p4],"bottom"],
#           [["spline3",p4,p5,p6],"bottom"],
#           [["spline3",p6,p7,p8],"bottom"],
#           [["line",p8,p9],"bottom"],
#           [["line",p9,p1],"left"]]
# [geo2.Append(c,bc=bc) for c,bc in curves]

# Creating the waveguide geometry.
# WG = MoveTo(0,0).Circle(0, 0, 1).Face() # Creating the rectangle of size WxD, and a circular scatterer at (x_sc,y_sc), of size b.
# PMLOUT = MoveTo(0,0).Circle(0, 0, 2).Face()

#Combining the 3 domains.
# geo = OCCGeometry(WG, dim=2) 

# geo.SetMaterial(1,"round")
# geo.SetMaterial(3,"PMLR")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))
mesh.Curve(3)
Draw(mesh)
print("Ok")
   
# Saving the mesh as Gmsh2 format.
meshname = "Mesh_saving_test.msh"
mesh.ngmesh.Export(meshname,"Gmsh2 Format") # Saving the mesh file. Not needed if you choose to interpolate to a regular grid later on.