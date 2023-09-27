# Preamble, loading necessary packages
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
from scipy.io import savemat
from scipy.interpolate import griddata
import numpy as np
import meshio



geo = SplineGeometry()
p1,p2,p3,p4 = [ geo.AppendPoint(x,y) for x,y in [(0,0), (1,0), (1,1), (0,1)] ]
p5,p6 =  [ geo.AppendPoint(x,y) for x,y in [(2,0), (2,1)] ]
geo.Append (["line", p1, p2], leftdomain=1, rightdomain=0)
geo.Append (["line", p2, p3], leftdomain=1, rightdomain=2)
geo.Append (["line", p3, p4], leftdomain=1, rightdomain=0)
geo.Append (["line", p4, p1], leftdomain=1, rightdomain=0)
geo.Append (["line", p2, p5], leftdomain=2, rightdomain=0)
geo.Append (["line", p5, p6], leftdomain=2, rightdomain=0)
geo.Append (["line", p6, p3], leftdomain=2, rightdomain=0)


mesh =Mesh(geo.GenerateMesh(maxh=0.2))
mesh.Curve(3)

mesh.SetPML(pml.Cartesian((0.0),(1,1),2j),2)
Draw(mesh)
