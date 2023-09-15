from ngsolve import *
#from ngsolve.webgui import Draw
from netgen.occ import *
from netgen.geom2d import SplineGeometry

'''
air = Rectangle(500, 200).Face()
air.edges.name = 'outer'
scatterer = Circle((350, 100),20).Face()
scatterer.edges.name = 'scat'
'''
#air = Rectangle(580,200).Circle(375,100,20).Reverse().Face()
air = MoveTo(-100,0).Rectangle(600,200).Face()
geo = OCCGeometry(air, dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=10))
mesh.Curve(3)
mesh.GetBoundaries
#Draw(mesh);

mesh.SetPML(pml.Cartesian((-100,0), (  0,200), 1j),"PMLLEFT")
mesh.SetPML(pml.Cartesian(( 500,0), (600,200), 1j),"PMLRIGHT")
#pc2 = pml.Cartesian((-40.,0.),(0, 200),5j)
#mesh.SetPML(pc2,"PMLL")
#pc3 = pml.Cartesian((540.,0.),(580, 200),1j)
#mesh.SetPML(pc3,"PMLR")

fes = H1(mesh, order=5, complex=True)
u, v = fes.TnT()

#Wavenumber & source
cspeed = 1500.
freq = 20.
x0=60.
y0=116.
omega = 2.*pi*freq / cspeed
ScaleFactor = 1e1

pulse = ScaleFactor *exp(-(omega**2)*((x-x0)*(x-x0) + (y-y0)*(y-y0)))

Draw(pulse, mesh, 'pulse')


a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx - omega*omega*u*v*dx
a.Assemble()

pulse = ScaleFactor *exp(-(omega**2)*((x-x0)*(x-x0) + (y-y0)*(y-y0)))
f = LinearForm(fes)
f += -pulse * v * dx
f.Assemble()
gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec

Draw(gfu)
