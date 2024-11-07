# wave equation with cartesian PML and DG
# Joachim Schöberl
# Preamble, loading necessary packages
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
from scipy.io import savemat
from scipy.interpolate import griddata
import numpy as np
import meshio
import os
import ngsolve.internal as ngsint
import netgen.gui



## Problem setup.
# Wavenumber & source.
c0 = 1500. # Constant wave speed.
f0 = 75.   # Reference frequency.
lambda_0 = c0/f0; # Reference wavelength. We use this to define all sizes with respect to this scale.

x_s = 50
y_s = 100
f = 73. # Frequency in which the source emits its pulse.
r = 30
alpha = np.log(10**6)/(r**2)
pulse = (np.sqrt(alpha/np.pi)*exp(-alpha*((x-x_s)*(x-x_s) + (y-y_s)*(y-y_s))) )
# Waveguide characteristics
D = 10*lambda_0 # Waveguide depth
W = 10*lambda_0 # Waveguide width. (originally infinite but we have to truncate for computational purposes).
PML_size = 4*lambda_0 # Length of the Perfectly Matched Layer (PML) that helps us truncate our computational domain.


geo = SplineGeometry()

pntx = [-PML_size,0,W,W+PML_size]
pnty = [0,D]
pmlpnts = []
for yi in pnty:
    for xi in pntx:
        pmlpnts.append (geo.AddPoint(xi,yi))
        print((xi,yi))
#inner rect
geo.Append (["line", pmlpnts[1], pmlpnts[2]], leftdomain=1, rightdomain=0)
geo.Append (["line", pmlpnts[2], pmlpnts[6]], leftdomain=1, rightdomain=2)
geo.Append (["line", pmlpnts[6], pmlpnts[5]], leftdomain=1, rightdomain=0)
geo.Append (["line", pmlpnts[5], pmlpnts[1]], leftdomain=1, rightdomain=2)

# right x-pml
geo.Append (["line", pmlpnts[2], pmlpnts[3]], leftdomain=2, rightdomain=0)
geo.Append (["line", pmlpnts[3], pmlpnts[7]], leftdomain=2, rightdomain=0)
geo.Append (["line", pmlpnts[7], pmlpnts[6]], leftdomain=2, rightdomain=0)
# left x-pml
geo.Append (["line", pmlpnts[5], pmlpnts[4]], leftdomain=2, rightdomain=0)
geo.Append (["line", pmlpnts[4], pmlpnts[0]], leftdomain=2, rightdomain=0)
geo.Append (["line", pmlpnts[0], pmlpnts[1]], leftdomain=2, rightdomain=0)

geo.SetMaterial(1,"inner")
geo.SetMaterial(2,"pmlx")


mesh = Mesh( geo.GenerateMesh(maxh=2.5))
Draw (mesh)

k = 4


fes_u = VectorL2(mesh, order=k)
fes_p = L2(mesh, order=k+1) 
fes = FESpace( [fes_p,fes_p,fes_u, fes_u] )
p,phat, u,uhat = fes.TrialFunction()
q,qhat, v,vhat = fes.TestFunction()

n = specialcf.normal(2) 

B  = BilinearForm(fes)
B += SymbolicBFI ( grad(p)*v )
B += SymbolicBFI ( -0.5 * (p - p.Other()) * (n*v), element_boundary = True)

Bt  = BilinearForm(fes)
Bt += SymbolicBFI ( -grad(q)*u )
Bt += SymbolicBFI ( 0.5 * (q-q.Other()) * (n*u), element_boundary = True) 

ux = u[0]
uy = u[1]
vx = v[0]
vy = v[1]
uxhat = uhat[0]
uyhat = uhat[1]
vxhat = vhat[0]
vyhat = vhat[1]

# damping matrices
sigma = 1
dampp = BilinearForm(fes)
dampp += SymbolicBFI ( p*q, definedon=mesh.Materials("pmlx"))

dampu = BilinearForm(fes)
dampu += SymbolicBFI ( ux*vx,  definedon=mesh.Materials("pmlx"))
dampu += SymbolicBFI ( -uy*vy + uy*vyhat-vy*uyhat + uyhat*vyhat,  definedon=mesh.Materials("pmlx"))

gfu = GridFunction(fes)
gfu.components[0].Set (pulse)
ngsint.viewoptions.drawoutline=0 # disable triangle outline when plotting.
Draw(gfu.components[0], mesh, "p",autoscale=False,min=-5e-3,max = 5e-3)
# Draw(gfu.components[1], mesh, "phat")
# Draw(gfu.components[2], mesh, "u")
# Draw(gfu.components[3], mesh, "uhat")

# For the time stepping
tau = 2e-6*c0*30
tend = 10*c0*5
t = 0

pdofs = fes.Range(0)
phatdofs = fes.Range(1)
udofs = fes.Range(2)
uhatdofs = fes.Range(3)

w = gfu.vec.CreateVector()
hv = gfu.vec.CreateVector()
i = 0
with TaskManager():
    while t < tend:
        print (t)
        Bt.Apply (gfu.vec, w)
        dampp.Apply (gfu.vec, hv)
        w.data -= sigma * hv
        fes_p.SolveM (vec=w[pdofs])
        fes_p.SolveM (vec=w[phatdofs])
        gfu.vec.data += tau * w

        B.Apply (gfu.vec, w)
        dampu.Apply (gfu.vec, hv)
        w.data -= sigma * hv
        fes_u.SolveM (vec=w[udofs])
        fes_u.SolveM (vec=w[uhatdofs])
        gfu.vec.data += tau * w

        t += tau
        i+=1
        Redraw()
        # Redraw(False)


