# Preamble, loading necessary packages
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
from scipy.io import savemat
from scipy.interpolate import griddata
import numpy as np
import meshio
from configs import *
import ngsolve.internal as ngsint
from tqdm import tqdm
ngsint.viewoptions.drawoutline=0 # disable triangle outline when plotting.


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
#
#           -----------------------------------------------------------------------------
#                                                             
#                                                             
#                                                        
#
#           ---------------------------
#                                       \
#                                        \                   
#                                         \
#                                          \
#                                           \
#                                            --------------------------------------------



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
geo.SetMaterial(2,"PMLL")
geo.SetMaterial(3,"PMLR")
mesh = Mesh(geo.GenerateMesh(maxh=5))
mesh.Curve(3)
#Draw(mesh)


# # Letting the mesh know which are the PML layers. Notice that we essentially redefine the WGBOX, add the absorbing parameter (2j).
# # We then tell the mesh which labels correspond to PMLs. In our case it's the PMLLEFT and PMLRIGHT faces we defined earlier.
mesh.SetPML(pml.Cartesian((0,0), (Wm,Dm), 2j),"PMLL|PMLR") 

# #Draw(mesh); Optional meshing drawing to see our work before proceeding.

# # Creating the finite element space based on the mesh we just created.
# # We define the order of the finite elements, and boundary conditions (leave blank for Neummann b.c.'s)
fes = H1(mesh, order=2, complex=True, dirichlet='top|bottom|scatterer') 

u, v = fes.TnT() # Creating Test and Trial functions u, v.


# # Source present in the waveguide.
# frq = 33. # Frequency in which the source emits its pulse.
Presp = np.zeros((Nr,Nr),dtype='complex')
Gsave = np.zeros((Nr, mesh.nv),dtype='complex')
# Gimw  = np.zeros((Nx,Ny),dtype='complex')
print("Generating incident field")
for n in tqdm(range(Nr)):
    y_s=  y_a[n] # Position of source in y-axis.
    pulse = sqrt(alpha/pi)*exp(-alpha*((x-x_s)*(x-x_s) + (y-y_s)*(y-y_s)))

    #Draw(pulse, mesh,'mesh') # Optional drawing to see what the source looks like.

    # Creating the weak form of the Helmholtz equation -Du - k^2 u = f
    a = BilinearForm(fes, symmetric=True) # Setting a as a bilinear form
    a += grad(u)*grad(v)*dx - k*k*u*v*dx
    a.Assemble()

    f = LinearForm(fes) # RHS is a linear form that contains the source.
    f += pulse * v * dx
    f.Assemble()

    # Create the grid function gfu that will contain our solution.
    gfu = GridFunction(fes, name="u")

    # Solve the linear system.
    gfu.vec.data = a.mat.Inverse() * f.vec

    # Draw the modulus of the complex solution on the mesh.
    # Draw(Norm(gfu),mesh,'mesh',)
    # scene = Draw(Norm(gfu),mesh,'mesh',autoscale=False,min=0,max=3,deformation=True)
#     Draw(Norm(gfu),mesh,"mesh", deformation=True, settings = {"camera" :{"transformations" :
#                                     [{ "type": "rotateY", "angle": -90},
#                                      { "type": "rotateX", "angle": 0}]}},
#      min=0, max=5, autoscale=False)
# print('Ok2')
    # Saving the mesh as Gmsh2 format.
    meshname = "mesh_inc_var_depth.msh"
    mesh.ngmesh.Export(meshname,"Gmsh2 Format") # Saving the mesh file. Not needed if you choose to interpolate to a regular grid later on.

# Saving the solution to a .mat file.
    sol_on_mesh = ConvertSolutiononMesh(mesh,gfu) # Only keeping parts of the solution that are on mesh points and not all DOFs.
    grid_x, grid_y = np.meshgrid(xg, yg) # Creating the regular grid we interpolate over.
    
# Making the meshpoints from ngmesh into a numpy array, in order to be able to use them on the griddata command.
# Kinda messy right now, quite possible doable in a better way.
    mesh_points =  np.zeros((mesh.nv,2))
    i = 0
    for p in mesh.ngmesh.Points():
        mesh_points[i,0] = p[0]
        mesh_points[i,1] = p[1]
        i = i + 1
    sol_on_array = griddata(mesh_points, sol_on_mesh, (x_s,y_a), method='cubic') # Interpolate the solution from the mesh points into the regular grid points.
    # sol_on_imw = griddata(mesh_points, sol_on_mesh, (grid_x,grid_y), method='cubic') # Interpolate the solution from the mesh points into the regular grid points.
    if not os.path.exists(os.path.join(os.path.dirname(__file__), 'data')): # If the data folder doesn't exist.
        os.mkdir(os.path.join(os.path.dirname(__file__), 'data'))       # Then this creates it.
    Presp[n,:] = sol_on_array
    Gsave[n,:] = sol_on_mesh 
    # Gsave[n,:,:] = sol_on_imw 
    
    
data_name = "data/inc_f"+ str(frq) + ".mat"
savemat(data_name,{"u":Presp}) # Save a mat file of the solution on a regular grid.

data_name = "data/green_f"+ str(frq) + ".mat"
savemat(data_name,{"u":Gsave}) # Save a mat file of the Green's function on a regular grid.

data_name = "data/mesh_points_inc.mat"
savemat(data_name,{"u":mesh_points})