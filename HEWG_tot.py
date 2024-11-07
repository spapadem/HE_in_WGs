# Preamble, loading necessary packages
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
from scipy.io import savemat
from scipy.interpolate import griddata
import numpy as np
import meshio
import ngsolve.internal as ngsint


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
#                                                               __
#                                                              /  \
#                                                              \__/
#
#           -----------------------------------------------------------------------------

## Problem setup.
# Wavenumber & source.
c0 = 1500. # Constant wave speed.
f0 = 75.   # Reference frequency.
lambda_0 = c0/f0; # Reference wavelength. We use this to define all sizes with respect to this scale.


# Waveguide characteristics
D = 10*lambda_0 # Waveguide depth
W = 25*lambda_0 # Waveguide width. (originally infinite but we have to truncate for computational purposes).
PML_size = 2*lambda_0 # Length of the Perfectly Matched Layer (PML) that helps us truncate our computational domain.

# Location and size of the scatterer.
x_sc = 19*lambda_0 # Location in x-axis
y_sc = 6*lambda_0 # Location in y-axis
b = 3*lambda_0 # Size of the scatterer (radius)

# Creating the waveguide geometry.
WG = MoveTo(0,0).Rectangle(W,D).Circle(x_sc, y_sc, b).Reverse().Face() # Creating the rectangle of size WxD, and a circular scatterer at (x_sc,y_sc), of size b.
WG.faces.name='WGBOX' # Naming the faces for the waveguide domain.
WG.edges.Min(Y).name='top' # Labeling the boundary at Y = 0.
WG.edges.Max(Y).name='bottom' # Labeling the boundary at Y = D.
WG.edges.Min(X).name='left' # Labeling the boundary at X = 0.
WG.edges.Max(X).name='right' # Labeling the boundary at X = W.

# Truncating the waveguide's infinite domain to a box using PMLs. 
# We will use one PML on the left (named PMLL) and one on the right (PMLR).
#           -----------------------------------------------------------------------------
#               |  P |                                  __             |  P  |
#               |  M |                                 /  \            |  M  |
#               |  L |                                 \__/            |  L  |
#               |  L |                                                 |  R  |
#           -----------------------------------------------------------------------------

# Defining PMLL.
PMLL = MoveTo(-PML_size,0).Rectangle(PML_size,D).Face() # The PMLL is a rectangle to the left side of the origin, spanning [-PML_size, 0] x [0 D]
PMLL.faces.name = 'PMLLEFT' # Naming the faces of the PMLL layer. We will need this for the SetPML command later.
PMLL.edges.Min(X).name = 'pml_left_b' # Labeling the boundary at X = -PML_size
PMLL.edges.Min(Y).name = 'top' # Labeling the boundary at Y = 0. We give it the same name as the WG domain, as we will need to apply the same b.c.'s.
PMLL.edges.Max(Y).name = 'bottom' # Labeling the boundary at Y = D. We give it the same name as the WG domain, as we will need to apply the same b.c.'s.

# Defining PMLR.
PMLR = MoveTo(W,0).Rectangle(PML_size,D).Face() #  The PMLL is a rectangle to the left side of the origin, spanning [W, W+PML_size] x [0 D]
PMLR.faces.name = 'PMLRIGHT' # Naming the faces of the PMLL layer. We will need this for the SetPML command later.
PMLR.edges.Max(X).name = 'pml_right_b' # Labeling the boundary at X = D + PML_size
PMLR.edges.Min(Y).name = 'top'# Labeling the boundary at Y = 0. We give it the same name as the WG domain, as we will need to apply the same b.c.'s.
PMLR.edges.Max(Y).name = 'bottom' # Labeling the boundary at Y = D. We give it the same name as the WG domain, as we will need to apply the same b.c.'s.

#Combining the 3 domains.
geo = OCCGeometry(Glue([PMLL,WG,PMLR]), dim=2) # Using the Glue command to combine PMLL, WG and PMLR. The order does not matter.

# Meshing
mesh_init = geo.GenerateMesh(maxh=5)
mesh = Mesh(mesh_init) # Creating a netgen mesh based on the geometry we just created.
mesh.Curve(3) # Defining mesh curve.

# Letting the mesh know which are the PML layers. Notice that we essentially redefine the WGBOX, add the absorbing parameter (2j).
# We then tell the mesh which labels correspond to PMLs. In our case it's the PMLLEFT and PMLRIGHT faces we defined earlier.
mesh.SetPML(pml.Cartesian((0,0), (W,D), 2j),"PMLLEFT|PMLRIGHT") 

#Draw(mesh); Optional meshing drawing to see our work before proceeding.

# Creating the finite element space based on the mesh we just created.
# We define the order of the finite elements, and boundary conditions (leave blank for Neummann b.c.'s)
fes = H1(mesh, order=5, complex=True, dirichlet='top|bottom|plm_left_b|pml_right_b') 

u, v = fes.TnT() # Creating Test and Trial functions u, v.


# Source present in the waveguide.
f = 173. # Frequency in which the source emits its pulse.
omega = 2.*pi*f / c0 # Angular frequency.
x_s=50. # Position of source in x-axis.
y_s=30. # Position of source in y-axis.
pulse = exp(-(omega**2)*((x-x_s)*(x-x_s) + (y-y_s)*(y-y_s)))


Draw(pulse, mesh,'mesh') # Optional drawing to see what the source looks like.

# Creating the weak form of the Helmholtz equation -Du - k^2 u = f
a = BilinearForm(fes, symmetric=True) # Setting a as a bilinear form
a += grad(u)*grad(v)*dx - omega*omega*u*v*dx
a.Assemble()

f = LinearForm(fes) # RHS is a linear form that contains the source.
f += -pulse * v * dx
f.Assemble()

# Create the grid function gfu that will contain our solution.
gfu = GridFunction(fes, name="u")

# Solve the linear system.
gfu.vec.data = a.mat.Inverse() * f.vec

# Draw the modulus of the complex solution on the mesh.
ngsint.viewoptions.drawoutline=0 # disable triangle outline when plotting.
Draw(Norm(gfu),mesh,'mesh')

# Saving the mesh as Gmsh2 format.
meshname = "Mesh_saving_test.msh"
mesh.ngmesh.Export(meshname,"Gmsh2 Format") # Saving the mesh file. Not needed if you choose to interpolate to a regular grid later on.

# Saving the solution to a .mat file.
sol_on_mesh = ConvertSolutiononMesh(mesh,gfu) # Only keeping parts of the solution that are on mesh points and not all DOFs.
Nx = 701 # Number of points in the regular grid, in the x-direction.
Ny = 201 # Number of points in the regular grid, in the y-direction.
grid_x, grid_y = np.meshgrid(np.linspace(-PML_size, W+PML_size, Nx),np.linspace(0, D, Ny)) # Creating the regular grid we interpolate over.

# Making the meshpoints from ngmesh into a numpy array, in order to be able to use them on the griddata command.
# Kinda messy right now, quite possible doable in a better way.
mesh_points =  np.zeros((mesh.nv,2))
i = 0
for p in mesh.ngmesh.Points():
     mesh_points[i,0] = p[0]
     mesh_points[i,1] = p[1]
     i = i + 1
sol_on_grid = griddata(mesh_points, sol_on_mesh, (grid_x, grid_y), method='cubic') # Interpolate the solution from the mesh points into the regular grid points.

if not os.path.exists(os.path.join(os.path.dirname(__file__), 'data')): # If the data folder doesn't exist.
        os.mkdir(os.path.join(os.path.dirname(__file__), 'data'))       # Then this creates it.
savemat("data/data_total.mat",{"u":sol_on_grid}) # Save a mat file of the solution on a regular grid.
