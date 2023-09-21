from ngsolve import *
#from ngsolve.webgui import Draw
from netgen.occ import *
from netgen.geom2d import SplineGeometry


D = 200;
W = 500;
PML_size = 20;


'''
air = Rectangle(500, 200).Face()
air.edges.name = 'outer'
scatterer = Circle((350, 100),20).Face()
scatterer.edges.name = 'scat'
'''


#air = MoveTo(-PML_size,0).Rectangle(W+PML_size,D).Circle(375,100,20).Reverse().Face()
air = MoveTo(0,0).Rectangle(W,D).Face()
air.faces.name='WGBOX'
air.edges.Min(Y).name='top'
air.edges.Max(Y).name='bottom'
air.edges.Min(X).name='left'
air.edges.Max(X).name='right'
geo = OCCGeometry(air, dim=2)
PMLL = MoveTo(-PML_size,0).Rectangle(PML_size,D).Face()
PMLL.faces.name = 'PMLLEFT'
PMLL.edges.Min(X).name = 'pml_left_b'
PMLL.edges.Min(Y).name = 'top'
PMLL.edges.Max(Y).name = 'bottom'
PMLR = MoveTo(W,0).Rectangle(PML_size,D).Face()
PMLR.faces.name = 'PMLRIGHT'
PMLR.edges.Max(X).name = 'pml_right_b'
PMLR.edges.Min(Y).name = 'top'
PMLR.edges.Max(Y).name = 'bottom'
geo = OCCGeometry(Glue([PMLL,air,PMLR]), dim=2)
mesh = Mesh(geo.GenerateMesh(maxh=5)) # this is a netgen mesh
mesh.Curve(3) #this is an NGsolve mesh
#Draw(mesh);

#mesh.SetPML(pml.Cartesian((-PML_size,0), (0,D), 100),"PMLLEFT")
#mesh.SetPML(pml.Cartesian((W,0), (W+PML_size,D), 100),"PMLRIGHT")
mesh.SetPML(pml.Cartesian((0,0), (W,D), 2j),"PMLLEFT|PMLRIGHT")

fes = H1(mesh, order=8, complex=True, dirichlet='top|bottom|plm_left_b|pml_right_b')
#fes = H1(mesh, order=5, complex=True)
u, v = fes.TnT()

#Wavenumber & source
cspeed = 1500.
freq = 70.
x0=125.
y0=116.
omega = 2.*pi*freq / cspeed
ScaleFactor = 1e1

pulse = ScaleFactor *exp(-(omega**2)*((x-x0)*(x-x0) + (y-y0)*(y-y0)))

Draw(pulse, mesh,'mesh')


a = BilinearForm(fes, symmetric=True)
a += grad(u)*grad(v)*dx - omega*omega*u*v*dx
a.Assemble()

pulse = ScaleFactor *exp(-(omega**2)*((x-x0)*(x-x0) + (y-y0)*(y-y0)))
f = LinearForm(fes)
f += -pulse * v * dx
f.Assemble()
gfu = GridFunction(fes, name="u")
gfu.vec.data = a.mat.Inverse() * f.vec

Draw(Norm(gfu),mesh,'mesh',draw_vol=False)
#Draw(m, objects=[lines,points, text], settings={"Objects": {"Edges": False, "Surface": False}})


vtk = VTKOutput(ma=mesh, coefs=[gfu.real], names = ["ureal"],filename="wgfem_real",subdivision=2)
vtk.Do()
vtk = VTKOutput(ma=mesh, coefs=[gfu.imag], names = ["uimag"],filename="wgfem_imag",subdivision=2)
vtk.Do()


# ngmesh.save(filename.vol)
# meshioi.read(filename)
# mesh.points
gfu.vec.data
