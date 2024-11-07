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
# geo = SplineGeometry()

# pnts =[(0,0), #1
#        (Wm,0), #2
#        (Wm,Dm),  #3
#        (Wm-4*lambda_0,Dm), #4
#        (Wm-7*lambda_0,Dm), #5
#        (0.5*(Wc+Wm-5*lambda_0), (Dc+Dm)/2), #6
#        (Wc,Dc), #7
#        (Wc-2*lambda_0,Dc), #8
#        (0,Dc)] #9
# p1,p2,p3,p4,p5,p6,p7,p8,p9 = [geo.AppendPoint(*pnt) for pnt in pnts]
# curves = [[["line",p1,p2],"top"],
#           [["line",p2,p3],"right"],
#           [["line",p3,p4],"bottom"],
#           [["spline3",p4,p5,p6],"bottom"],
#           [["spline3",p6,p7,p8],"bottom"],
#           [["line",p8,p9],"bottom"],
#           [["line",p9,p1],"left"]]
# [geo.Append(c,bc=bc) for c,bc in curves]

# Creating the waveguide geometry.
WG = MoveTo(0,0).Circle(0, 0, 1).Face() # Creating the rectangle of size WxD, and a circular scatterer at (x_sc,y_sc), of size b.
# PMLOUT = MoveTo(0,0).Circle(0, 0, 2).Face()

#Combining the 3 domains.
geo = OCCGeometry(WG, dim=2) 

# geo.SetMaterial(2,"PMLL")
# geo.SetMaterial(3,"PMLR")
mesh = Mesh(geo.GenerateMesh(maxh=0.1))
mesh.Curve(3)
Draw(mesh)
print("Ok")
   
# Saving the mesh as Gmsh2 format.
meshname = "Mesh_saving_test.msh"
mesh.ngmesh.Export(meshname,"Gmsh2 Format") # Saving the mesh file. Not needed if you choose to interpolate to a regular grid later on.