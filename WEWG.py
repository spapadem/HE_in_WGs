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

# This function saves the solution that typically lives on all the DOFs, just on the mesh points. 
# Useful for visualization purposes when using higher order elements, as the exported mesh seems
# to not include the DOFs as well. We will use this function to save the solution on just the mesh points
# and then interpolate on a regular grid.
def ConvertSolutiononMesh(mesh,gfu):
    soln = np.zeros(mesh.nv,dtype='complex')
   
    i = 0
    for p in mesh.ngmesh.Points():
        soln[i] = gfu(p[0],p[1])
        i = i + 1
    return soln

# We work on an infinite waveguide, which means it is bounded on top and bottom and the solution is outgoing on the left and right.
# Specifically, we will consider an infinite strip waveguide, with flat top and bottom boundaries. A waveguide with a varying bottom
# will be implemented in a future project.
#
#           -----------------------------------------------------------------------------
#
#
#
#
#           -----------------------------------------------------------------------------

# Wavenumber & source.
c0 = 1500. # Constant wave speed.
f0 = 75.   # Reference frequency.
lambda_0 = c0/f0; # Reference wavelength. We use this to define all sizes with respect to this scale.


# Waveguide characteristics
Dc = 7*lambda_0 # Waveguide constant depth
Dm = 20*lambda_0 # Waveguide max depth
Wc = 10*lambda_0 # Width with constant depth
Wm = 25*lambda_0 # Waveguide max width. (originally infinite but we have to truncate for computational purposes).
PML_size = 4*lambda_0 # Length of the Perfectly Matched Layer (PML) that helps us truncate our computational domain.

# Location and size of the scatterer.
x_sc = 20*lambda_0 # Location in x-axis
y_sc = 3.5*lambda_0 # Location in y-axis
b = 2*lambda_0 # Size of the scatterer (radius)

# Creating the waveguide geometry.
geo = SplineGeometry()

pnts =[(0,0), #1
       (Wm,0), #2
       (Wm,Dm),  #3
       (Wm-4*lambda_0,Dm), #4
       (Wm-7*lambda_0,Dm), #5
       (0.5*(Wc+Wm-5*lambda_0), (Dc+Dm)/2), #6
       (Wc,Dc), #7
       (Wc-2*lambda_0,Dc), #8
       (0,Dc)] #9
p1,p2,p3,p4,p5,p6,p7,p8,p9 = [geo.AppendPoint(*pnt) for pnt in pnts]
curves = [[["line",p1,p2],"top"],
          [["line",p2,p3],"right"],
          [["line",p3,p4],"bottom"],
          [["spline3",p4,p5,p6],"bottom"],
          [["spline3",p6,p7,p8],"bottom"],
          [["line",p8,p9],"bottom"],
          [["line",p9,p1],"left"]]
[geo.Append(c,bc=bc) for c,bc in curves]

geo.AddRectangle((-PML_size,0),(0,Dc),leftdomain=2,bc="PMLL") # Add left PML rectangle.
geo.AddRectangle((Wm,0),(Wm+PML_size,Dm),leftdomain=3,bc="PMLR") # Add right PML rectangle.
geo.AddCircle((x_sc,y_sc),b,leftdomain=0,rightdomain=1,bc="scatterer") # Add scatterer in the domain.
geo.SetMaterial(2,"PMLL")
geo.SetMaterial(3,"PMLR")
mesh = Mesh(geo.GenerateMesh(maxh=5))
mesh.Curve(3)
#Draw(mesh)
print("Ok")

# # Letting the mesh know which are the PML layers. Notice that we essentially redefine the WGBOX, add the absorbing parameter (2j).
# # We then tell the mesh which labels correspond to PMLs. In our case it's the PMLLEFT and PMLRIGHT faces we defined earlier.

mesh.SetPML(pml.Cartesian((0,0), (Wm,Dm), 2j),"PMLL|PMLR") 

# #Draw(mesh); Optional meshing drawing to see our work before proceeding.

# # Creating the finite element space based on the mesh we just created.
# # We define the order of the finite elements, and boundary conditions (leave blank for Neummann b.c.'s)
fes = H1(mesh, order=3, complex=False, dirichlet='top|bottom|scatterer') 

u, v = fes.TnT() # Creating Test and Trial functions u, v.

# Source setup
x_s = 100 # Source location in x.
y_s = 50 # Source location in y.
r = 20 # Radius of source
alpha = np.log(10**6)/(r**2) # Normalization parameter
pulse = np.sqrt(alpha/np.pi)*exp(-alpha*((x-x_s)*(x-x_s) + (y-y_s)*(y-y_s))) # Gaussian pulse as an initial condition
# pulse = x*y*(2-x)*(3-y)

# Creating the weak form of the wave equation.
M = BilinearForm(fes, symmetric=True) # Setting Mass matrix.
M += -u*v*dx
M.Assemble()

#Inverting mass matrix for time-stepping
Mmat = M.mat.CreateMatrix()
Mmat.AsVector().data = M.mat.AsVector()
Minv = Mmat.Inverse(freedofs=fes.FreeDofs())


K = BilinearForm(fes, symmetric=True) # Setting Stiffness matrix.
K += -grad(u)*grad(v)*dx
K.Assemble()
Kmat = K.mat.CreateMatrix()
Kmat.AsVector().data = K.mat.AsVector()

# Create the grid function that will contain our solution.
u_sol = GridFunction(fes, name="u_sol")

u_old = GridFunction(fes, name="u_old")
u_old.Set(pulse)

u_older = GridFunction(fes, name="u_older")
u_older.Set(pulse)

# Plotting the solution
ngsint.viewoptions.drawoutline=0 # disable triangle outline when plotting.
Draw(u_old,mesh,'mesh',autoscale=False,min=-5e-3,max = 5e-3)

# Time-stepping
t = 0
T = 10
Dt = 0.000125
with TaskManager():
    while t <= T:
        print ("t=", t, end="\r")
        u_sol.vec.data = -Dt**2*c0**2*Minv*(K.mat*u_old.vec) + 2*u_old.vec - u_older.vec
        t = t + Dt
        u_older.vec.data = u_old.vec.data
        u_old.vec.data = u_sol.vec.data
        Redraw()
