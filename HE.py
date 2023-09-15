# solve the Poisson equation -Delta u = f
# with Dirichlet boundary condition u = 0

from ngsolve import *
from netgen.occ import *
#import numpy as np
#import matplotlib.pyplot as plt

# generate a triangular mesh of mesh-size 0.2

shape = MoveTo(-40,0).Rectangle(580,200).Circle(375,100,20).Reverse().Face()
shape.edges.name="wall"
shape.edges.Min(X).name="inlet"
shape.edges.Max(X).name="outlet"

#PML_left = MoveTo(-40,0).Rectangle(40,200).Face()
#PML_left.edges.name="pml_left"

#PML_right = MoveTo(500,0).Rectangle(40,200).Face()
#PML_right.edges.name="pml_right"



geo = OCCGeometry(shape, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=8))
mesh.Curve(3)
Draw (mesh);

#mesh.SetPML(pml.Cartesian((-20,0), (  0,200), 1j))
#mesh.SetPML(pml.Cartesian((500,0), (540,200), 1j))



#Wavenumber & source
cspeed = 1500.
freq = 20.
x0=20.
y0=36.
omega = 2.*pi*freq / cspeed
ScaleFactor = 1e1


fes = H1(mesh, order=2, complex=True)
u, v = fes.TnT()

a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx - omega*omega*u*v*dx
a.Assemble()

pulse = ScaleFactor *exp(-(omega**2)*((x-x0)*(x-x0) + (y-y0)*(y-y0)))
f = LinearForm(fes)
f += -pulse * v * dx
f.Assemble()
gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec


