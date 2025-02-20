# %% [markdown]
# # Spherical indenter.
# - Code for 3d large deformation viscoelasticity of elastomeric materials.
# - A cube under pressure on a part of its face.
# - Uses quadrature representation of internal variable: viscous deformation tensor Cv.
# 

# %% [markdown]
# ### Units
# - Basic:
#     - Length: mm
#     - Mass: kg  
#     - Time:  s
#     - Mass density: kg/mm^3
# 
# - Derived:
#     - Force: mN
#     - Stress: kPa 
#     - Energy: uJ
# 
# ### Software:
# - Dolfinx v0.8.0
#
# In the collection "Example Codes for Coupled Theories in Solid Mechanics,"
#
# By Eric M. Stewart, Shawn A. Chester, and Lallit Anand.
#
# https://solidmechanicscoupledtheories.github.io/

# %% [markdown]
# # Import modules

# %%
# Import FEnicSx/dolfinx
import dolfinx

# For numerical arrays
import numpy as np

# For MPI-based parallelization
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# PETSc solvers
from petsc4py import PETSc

# specific functions from dolfinx modules
from dolfinx import fem, mesh, io, plot, log
from dolfinx.fem import (Constant, dirichletbc, Function, functionspace, Expression )
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.io import VTXWriter, XDMFFile


# specific functions from ufl modules
import ufl
from ufl import (TestFunctions, TrialFunction, Identity, grad, det, div, dev, inv, tr, sqrt, conditional ,\
                 gt, dx, inner, derivative, dot, ln, split, exp, eq, cos, acos, ge, le)

# basix finite elements
import basix
from basix.ufl import element, mixed_element, quadrature_element

# Matplotlib for plotting
import matplotlib.pyplot as plt
plt.close('all')

# For timing the code
from datetime import datetime

# this forces the program to still print (but only from one CPU) 
# when run in parallel.
def mprint(*argv):
    if rank==0:
        out=""
        for arg in argv:
            out = out+ str(argv)
        print(out, flush=True)


# Set level of detail for log messages (integer)
# Guide:
# CRITICAL  = 50, // errors that may lead to data corruption
# ERROR     = 40, // things that HAVE gone wrong
# WARNING   = 30, // things that MAY go wrong later
# INFO      = 20, // information of general interest (includes solver info)
# PROGRESS  = 16, // what's happening (broadly)
# TRACE     = 13, // what's happening (in detail)
# DBG       = 10  // sundry
#
log.set_log_level(log.LogLevel.WARNING)

# %% [markdown]
# # Define geometry

# %%

# A 3-D cube
length = 50.0 # mm

# import the mesh and facet tags
with XDMFFile(MPI.COMM_WORLD,"meshes/sphere_indent.xdmf",'r') as infile:
    domain = infile.read_mesh(name="Grid",xpath="/Xdmf/Domain")
domain.topology.create_connectivity(domain.topology.dim, domain.topology.dim-1)

with XDMFFile(MPI.COMM_WORLD, "meshes/facet_sphere_indent.xdmf", "r") as xdmf:
    facet_tags = xdmf.read_meshtags(domain, name="Grid")

# This says "spatial coordinates" but is really the referential coordinates,
# since the mesh does not convect in FEniCS.
x = ufl.SpatialCoordinate(domain) 

# radius of indenter
r_ind = 10 # mm

# %% [markdown]
# **Print out the unique facet index numbers**

# # %%
# top_imap = domain.topology.index_map(2)      # index map of 2D entities in domain (facets)
# values = np.zeros(top_imap.size_global)      # an array of zeros of the same size as number of 2D entities
# values[facet_tags.indices]=facet_tags.values # populating the array with facet tag index numbers
# print(np.unique(facet_tags.values))          # printing the unique indices

# Surface numbering:
# boundaries = [(13,xBot),(14,xTop),(15,yBot),(16,yTop),(17,zBot),(18,zTop)]

# %% [markdown]
# **Visualize reference configuration and boundary facets**

# %%
# import pyvista
# pyvista.set_jupyter_backend('html')
# from dolfinx.plot import vtk_mesh
# pyvista.start_xvfb()

# # initialize a plotter
# plotter = pyvista.Plotter()

# # Add the mesh -- I make the 3D mesh opaque, so that 2D surfaces stand out.
# topology, cell_types, geometry = plot.vtk_mesh(domain, domain.topology.dim)
# grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
# plotter.add_mesh(grid, show_edges=True, opacity=0.25) 

# # Add colored 2D surfaces for the named surfaces
# xBot_surf = pyvista.UnstructuredGrid(*vtk_mesh(domain, domain.topology.dim-1,facet_tags.indices[facet_tags.values==13]) )
# yBot_surf = pyvista.UnstructuredGrid(*vtk_mesh(domain, domain.topology.dim-1,facet_tags.indices[facet_tags.values==15]) )
# zBot_surf = pyvista.UnstructuredGrid(*vtk_mesh(domain, domain.topology.dim-1,facet_tags.indices[facet_tags.values==17]) )
# #
# actor = plotter.add_mesh(xBot_surf, show_edges=True,color="red") # xBot face is blue
# actor2 = plotter.add_mesh(yBot_surf, show_edges=True,color="green") # yBot is red
# actor3 = plotter.add_mesh(zBot_surf, show_edges=True,color="blue") # zBot is green

# labels = dict(zlabel='Z', xlabel='X', ylabel='Y')
# plotter.add_axes(**labels)
    
# plotter.screenshot("results/cube_footing_mesh.png")

# from IPython.display import Image
# Image(filename='results/cube_footing_mesh.png') 

# %% [markdown]
# ### Un-comment this cell to see an interactive visualization of the mesh

# %%
# plotter.show()

# %% [markdown]
# ## Define boundary and volume integration measure

# %%
# Define the boundary integration measure "ds" using the facet tags,
# also specify the number of surface quadrature points.
ds = ufl.Measure('ds', domain=domain, subdomain_data=facet_tags, metadata={'quadrature_degree':2})

# Define the volume integration measure "dx" 
# also specify the number of volume quadrature points.
dx = ufl.Measure('dx', domain=domain, metadata={'quadrature_degree': 2, 'quadrature_rule':'default'})

#  Define facet normal
n = ufl.FacetNormal(domain)

# %% [markdown]
# # Material parameters

# %%
Geq_0   = Constant(domain, 400.0)      # Shear modulus, kPa
Kbulk = Constant(domain, PETSc.ScalarType(1.0e3*Geq_0))            # Bulk modulus, kPa
lambdaL = Constant(domain, 10.0)       # Arruda-Boyce locking stretch

# Viscoelasticity parameters from Linder et al. (2011) for NBR
#
Gneq_1  = Constant(domain, 535.7)    #  Non-equilibrium shear modulus, kPa
tau_1   = Constant(domain, 1.0)    #  relaxation time, s
#
Gneq_2  = Constant(domain, 76.20)    #  Non-equilibrium shear modulus, kPa
tau_2   = Constant(domain, 10.0)     #  relaxation time, s
# 
Gneq_3  = Constant(domain, 120.5)    #  Non-equilinrium shear modulus, kPa
tau_3   = Constant(domain, 100.0)    #  relaxation time, s
#
Gneq_4  = Constant(domain, 21.3)    #  Non-equilinrium shear modulus, kPa
tau_4   = Constant(domain, 1000.0)    #  relaxation time, s
#
Gneq_5  = Constant(domain, 22.9)    #  Non-equilinrium shear modulus, kPa
tau_5   = Constant(domain, 10000.0)    #  relaxation time, s

rho = Constant(domain, 1.3e-6) # Mass density, kg/mm^3

# alpha-method parameters
alpha   = Constant(domain, PETSc.ScalarType(0.0))
gamma   = Constant(domain, PETSc.ScalarType(0.5+alpha))
beta    = Constant(domain, PETSc.ScalarType(0.25*(gamma+0.5)**2))

# %% [markdown]
# # Simulation time-control related params

# %%
# Simulation time control-related params
t     = 0.0          # start time
tRamp = 60             # time for application of pressure 
Ttot  = 460           # total simulation time 
ind_max = 10    # depth of indentation, mm
numSteps = 50 #6*24     # number of steps for application of pressure ( larger number of steps can be used for smoother movies )
numSteps2 = 10 #5*24     # number of steps for relaxation after pressure application
dt = tRamp/numSteps  #  fixed step size
dt2 = (Ttot-tRamp)/numSteps2
dk = Constant(domain, PETSc.ScalarType(dt))

# Function to linearly ramp up indenter displacement
def dispRamp(t):
    
    if t < tRamp:
        disp = ind_max*t/tRamp 
    else:
        disp = ind_max
        
    return disp

# %% [markdown]
# # Function spaces

# %%

U2 = element("Lagrange", domain.basix_cell(), 2, shape=(3,))  # For displacement
P1 = element("Lagrange", domain.basix_cell(), 1)  # For pressure
# P0 = element("Lagrange", domain.basix_cell(), 1)
# T0 = element("Lagrange", domain.basix_cell(), 1, shape=(3,3))
P0 = quadrature_element(domain.basix_cell(), degree=2, scheme="default") 
# Note: it seems that for the current version of dolfinx, 
# only degree=2 quadrature elements actually function properly 
# in e.g. visualization interpolations and problem solution.
T0 = basix.ufl.blocked_element(P0, shape=(3, 3)) # for Cv
#
TH = mixed_element([U2, P1])     # Taylor-Hood style mixed element
ME = functionspace(domain, TH)    # Total space for all DOFs
#
V1 = functionspace(domain, P0) # Scalar function space.
V2 = functionspace(domain, U2) # Vector function space
V3 = functionspace(domain, T0) # Tensor function space
#
# Define actual functions with the required DOFs
w        = Function(ME)
u, p = split(w)  

# A copy of functions to store values in the previous step
w_old             = Function(ME)
u_old,  p_old = split(w_old)   

# Define test functions        
u_test, p_test = TestFunctions(ME)    

# Define trial functions needed for automatic differentiation
dw     = TrialFunction(ME)  
du, dp = split(dw)


# Define a tensor-valued function for Cv.
Cv_1_old = Function(V3) 
Cv_2_old = Function(V3) 
Cv_3_old = Function(V3)
Cv_4_old = Function(V3) 
Cv_5_old = Function(V3) 

# Functions for storing the velocity and acceleration at prev. step
v_old = Function(V2)
a_old = Function(V2)

# Initial conditions: 

# A function for constructing the identity matrix.
#
# To use the interpolate() feature, this must be defined as a 
# function of x.
def identity(x):
    values = np.zeros((domain.geometry.dim*domain.geometry.dim,
                      x.shape[1]), dtype=np.float64)
    values[0] = 1
    values[4] = 1
    values[8] = 1
    return values

# interpolate the identity onto the tensor-valued Cv function.
Cv_1_old.interpolate(identity)  
Cv_2_old.interpolate(identity) 
Cv_3_old.interpolate(identity)   
Cv_4_old.interpolate(identity) 
Cv_5_old.interpolate(identity)   

# %% [markdown]
# # Subroutines for kinematics and constitutive equations

# %%
#------------------------------------------------------------- 
# Utility subroutines
#-------------------------------------------------------------
 
# Subroutine for a "safer" sqrt() function which avoids a divide by zero 
# when differentiated. 
def safe_sqrt(x):
    return sqrt(x + 1.0e-16)

# Deformation gradient 
def F_calc(u):
    Id = Identity(3) 
    F = Id + grad(u) 
    return F

#------------------------------------------------------------- 
# Subroutines for computing the viscous flow update
#-------------------------------------------------------------

# subroutine for the distortional part / unimodular part of a tensor A
def dist_part(A):

    Abar = A / (det(A)**(1.0/3.0))

    return Abar

# Subroutine for computing the viscous stretch Cv at the end of the step.
def Cv_update(u, Cv_old, tau_r):
    
   F = F_calc(u)
   
   J = det(F)
   
   C = F.T*F
   
   Cv_new = dist_part( Cv_old + ( dk / tau_r ) * J**(-2./3.) * C ) 
    
   return Cv_new

#------------------------------------------------------------- 
# Subroutines for calculating the equilibrium Cauchy stress
#-------------------------------------------------------------

# Subrountine for computing the effective stretch
def lambdaBar_calc(u):
    
    F = F_calc(u)
    
    J = det(F)
    
    Fbar = J**(-1/3)*F
    
    Cbar = Fbar.T*Fbar
    
    I1 = tr(Cbar)
    
    lambdaBar = safe_sqrt(I1/3.0)
    
    return lambdaBar

# Subroutine for computing the zeta-function in the Arruda-Boyce model.
def zeta_calc(u):
    
    lambdaBar = lambdaBar_calc(u)
    
    # Use Pade approximation of Langevin inverse (A. Cohen, 1991)
    # This is sixth-order accurate.
    z    = lambdaBar/lambdaL
    
    z    = conditional(gt(z,0.99), 0.99, z) # Keep from blowing up
    
    beta = z*(3.0 - z**2.0)/(1.0 - z**2.0)
    
    zeta = (lambdaL/(3*lambdaBar))*beta
    
    return zeta

# Generalized shear modulus for Arruda-Boyce model
def Geq_AB_calc(u):
    
    zeta = zeta_calc(u)
    
    Geq_AB  = Geq_0 * zeta
    
    return Geq_AB

# Subroutine for calculating the  equilibrium Cauchy stress
def T_eq_calc(u,p):
    
    F   = F_calc(u)
    
    J = det(F)
    
    Fbar = J**(-1./3.)*F
    
    Bbar = Fbar*Fbar.T
    
    Geq  = Geq_AB_calc(u)
    
    T_eq = (1/J)* Geq * dev(Bbar) - p * Identity(3)
    
    return T_eq


#------------------------------------------------------------- 
# Subroutine for calculating the total Piola stress
#-------------------------------------------------------------
   
# Subroutine for the non-equilibrium Cauchy stress.
def T_neq_calc(u, Cv, Gneq):
        
    F  = F_calc(u)
    
    J = det(F)
    
    C = F.T*F
    
    T_neq = J**(-5./3.) * Gneq * (F * inv(Cv) * F.T - (1./3.) * inner(C, inv(Cv)) * Identity(3) ) 
    
    return T_neq
    
    
# Subroutine for the total Piola stress.
def  Piola_calc(u, p, Cv_1, Cv_2, Cv_3, Cv_4, Cv_5, Gneq_1, Gneq_2, Gneq_3, Gneq_4, Gneq_5):
    
    F  = F_calc(u)
    
    J = det(F)
    
    T_eq = T_eq_calc(u,p)
    
    T_neq_1 = T_neq_calc(u, Cv_1, Gneq_1)
    
    T_neq_2 = T_neq_calc(u, Cv_2, Gneq_2)
        
    T_neq_3 = T_neq_calc(u, Cv_3, Gneq_3)
    
    T_neq_4 = T_neq_calc(u, Cv_4, Gneq_4)
        
    T_neq_5 = T_neq_calc(u, Cv_5, Gneq_5)
    
    T = T_eq + T_neq_1 + T_neq_2 + T_neq_3 + T_neq_4 + T_neq_5
    
    Piola = J * T * inv(F.T)
    
    return Piola 

#---------------------------------------------------------------------
# Subroutine for updating  acceleration using the Newmark beta method:
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
#---------------------------------------------------------------------
def update_a(u, u_old, v_old, a_old):
    return (u-u_old-dk*v_old)/beta/dk**2 - (1-2*beta)/2/beta*a_old

#---------------------------------------------------------------------
# Subroutine for updating  velocity using the Newmark beta method
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
#---------------------------------------------------------------------
def update_v(a, u_old, v_old, a_old):
    return v_old + dk*((1-gamma)*a_old + gamma*a)

#---------------------------------------------------------------------
# alpha-method averaging function
#---------------------------------------------------------------------
def avg(x_old, x_new, alpha):
    return alpha*x_old + (1-alpha)*x_new

# Macaulay bracket function
def ppos(x):
    return (x+abs(x))/2.


# %% [markdown]
# # Evaluate kinematics and constitutive relations

# %%
# Get acceleration and velocity at end of step
a_new = update_a(u, u_old, v_old, a_old)
v_new = update_v(a_new, u_old, v_old, a_old)

# get avg (u,p) fields for generalized-alpha method
u_avg  = avg(u_old, u, alpha)
p_avg  = avg(p_old, p, alpha)

# Kinematical quantities
F  = F_calc(u_avg)
J  = det(F)
lambdaBar = lambdaBar_calc(u_avg)

# update the Cv tensors
Cv_1 = Cv_update(u_avg, Cv_1_old, tau_1)
Cv_2 = Cv_update(u_avg, Cv_2_old, tau_2)
Cv_3 = Cv_update(u_avg, Cv_3_old, tau_3)
Cv_4 = Cv_update(u_avg, Cv_4_old, tau_4)
Cv_5 = Cv_update(u_avg, Cv_5_old, tau_5)

#  Evaulate the total Piola stress
Piola = Piola_calc(u_avg, p_avg, Cv_1, Cv_2, Cv_3, Cv_4, Cv_5, Gneq_1, Gneq_2, Gneq_3, Gneq_4, Gneq_5)

# spatial coordinates, for use in contact penalty
u_var = ufl.variable(u)
x_sp  = ufl.variable(x + u_var)

# %% [markdown]
# # Weak forms

# %%

# initial position of center of indenter 
x3_ind = Constant(domain, PETSc.ScalarType(length + r_ind - dispRamp(0)))
x_ind = ufl.as_vector([0, 0, x3_ind])
#
# the penalty energy for contact enforcement
k_pen = Constant(domain, PETSc.ScalarType(1.0e2))
contact_energy = 0.5 * k_pen * ( ppos(r_ind**2 - (x_sp - x_ind)**2 )**2 ) 
contact_traction = ufl.diff(contact_energy, u_var)
#
# The weak form for the equilibrium equation + contact penalty enforcement.
Res_1  =  inner( Piola, grad(u_test))*dx \
          + dot(contact_traction, u_test)*ds
          #+ inner(rho * a_new, u_test)*dx 
          # (we exclude inertial effects for now)
              
# The auxiliary equation for the pressure
#
Res_2 = inner((J-1) + p_avg/Kbulk, p_test)*dx

# The weak form for the viscous strain increment updates

# The total residual
Res = Res_1 + Res_2 

# Automatic differentiation tangent:
a = derivative(Res, w, dw)

# %% [markdown]
# # Set-up output files

# %%
# results file name
results_name = "3D_cube_sphere_indent_MPI"

# Function space for projection of results
P1 = element("Lagrange", domain.basix_cell(), 1)
VV1 = fem.functionspace(domain, P1) # linear scalar function space
#
U1 = element("Lagrange", domain.basix_cell(), 1, shape=(3,)) 
VV2 = fem.functionspace(domain, U1) # linear Vector function space
#
T1 = element("Lagrange", domain.basix_cell(), 1, shape=(3,3)) 
VV3 = fem.functionspace(domain, T1) # linear tensor function space

# For visualization purposes, we need to re-project the stress tensor onto a linear function space before 
# we write it (and its components and the von Mises stress, etc) to the VTX file. 
#
# This is because the stress is a complicated "mixed" function of the (quadratic Lagrangian) displacements
# and the (quadrature representation) plastic strain tensor and scalar equivalent plastic strain. 
#
# First, define a function for setting up this kind of projection problem for visualization purposes:
def setup_projection(u, V):

    trial = ufl.TrialFunction(V)
    test  = ufl.TestFunction(V)   

    a = ufl.inner(trial, test)*dx
    L = ufl.inner(u, test)*dx

    projection_problem = dolfinx.fem.petsc.LinearProblem(a, L, [], \
        petsc_options={"ksp_type": "cg", "ksp_rtol": 1e-16, "ksp_atol": 1e-16, "ksp_max_it": 1000})
    
    return projection_problem

# Create a linear problem for projecting the stress tensor onto the linear tensor function space VV3.
#
tensor_projection_problem = setup_projection(Piola, VV3)
Piola_temp = tensor_projection_problem.solve()

# %%
# primary fields to write to output file
u_vis      = Function(VV2, name="disp")
p_vis      = Function(VV1, name="p")

# %%

# Mises stress
T     = Piola_temp*F.T/J
T0    = T - (1/3)*tr(T)*Identity(3)
Mises = sqrt((3/2)*inner(T0, T0))
Mises_vis= Function(VV1,name="Mises")
Mises_expr = Expression(Mises,VV1.element.interpolation_points())

# Cauchy stress components
T11 = Function(VV1)
T11.name = "T11"
T11_expr = Expression(T[0,0],VV1.element.interpolation_points())

T22 = Function(VV1)
T22.name = "T22"
T22_expr = Expression(T[1,1],VV1.element.interpolation_points())

T33 = Function(VV1)
T33.name = "T33"
T33_expr = Expression(T[2,2],VV1.element.interpolation_points())


# %%
# Effective stretch
lambdaBar_vis      = Function(VV1)
lambdaBar_vis.name = "LambdaBar"
lambdaBar_expr     = Expression(lambdaBar, VV1.element.interpolation_points())

# Volumetric deformation
J_vis      = Function(VV1)
J_vis.name = "J"
J_expr     = Expression(J, VV1.element.interpolation_points())

# %%
# set up the output VTX files.
file_results = VTXWriter(
    MPI.COMM_WORLD,
    "results/" + results_name + ".bp",
    [  # put the functions here you wish to write to output
        u_vis, p_vis, # DOF outputs
        Mises_vis, T11, T22, T33, # stress outputs
        lambdaBar_vis, J_vis, # Kinematical outputs
    ],
    engine="BP4",
)

def writeResults(t):
    
    # Update the output fields before writing to VTX.
    #
    u_vis.interpolate(w.sub(0))
    p_vis.interpolate(w.sub(1))
    #
    # re-project to smooth visualization of quadrature functions
    # before interpolating.
    Piola_temp = tensor_projection_problem.solve()
    Mises_vis.interpolate(Mises_expr)
    T11.interpolate(T11_expr)
    T22.interpolate(T22_expr)
    T33.interpolate(T33_expr)
    #
    lambdaBar_vis.interpolate(lambdaBar_expr)
    J_vis.interpolate(J_expr)
       
    # Finally, write output fields to VTX.
    #
    file_results.write(t) 

# %% [markdown]
# # Infrastructure for pulling out time history data (force, displacement, etc.)

# %%
# infrastructure for evaluating functions at a certain point efficiently
pointForEval = np.array([0, 0, length])
bb_tree = dolfinx.geometry.bb_tree(domain,domain.topology.dim)
cell_candidates = dolfinx.geometry.compute_collisions_points(bb_tree, pointForEval)
colliding_cells = dolfinx.geometry.compute_colliding_cells(domain, cell_candidates, pointForEval).array

# computing the reaction force using the stress field
traction = dot(Piola_temp, n)
Force    = dot(traction, n)*ds(18)
rxnForce = fem.form(Force) 

# %% [markdown]
# ## Boundary condtions

# %%
# Surface numbering:
# boundaries = [(13,xBot),(14,xTop),(15,yBot),(16,yTop),(17,zBot),(18,zTop)]

# Find the specific DOFs which will be constrained.
xBtm_u1_dofs = fem.locate_dofs_topological(ME.sub(0).sub(0), facet_tags.dim, facet_tags.find(13))
yBtm_u2_dofs = fem.locate_dofs_topological(ME.sub(0).sub(1), facet_tags.dim, facet_tags.find(15))
zBtm_u3_dofs = fem.locate_dofs_topological(ME.sub(0).sub(2), facet_tags.dim, facet_tags.find(17))
# #
# load_u1_dofs = fem.locate_dofs_topological(ME.sub(0).sub(0), facet_tags.dim, facet_tags.find(18))
# load_u2_dofs = fem.locate_dofs_topological(ME.sub(0).sub(1), facet_tags.dim, facet_tags.find(18))

# building Dirichlet BCs
bcs_1 = dirichletbc(0.0, xBtm_u1_dofs, ME.sub(0).sub(0))  # u1 fix - xBtm
bcs_2 = dirichletbc(0.0, yBtm_u2_dofs, ME.sub(0).sub(1))  # u2 fix - yBtm
bcs_3 = dirichletbc(0.0, zBtm_u3_dofs, ME.sub(0).sub(2))  # u3 fix - zBtm
#
# bcs_4 = dirichletbc(0.0, load_u1_dofs, ME.sub(0).sub(0))  # u1 fix - xBtm
# bcs_5 = dirichletbc(0.0, load_u2_dofs, ME.sub(0).sub(1))  # u2 fix - yBtm

bcs = [bcs_1, bcs_2, bcs_3] #, bcs_4, bcs_5]

# %% [markdown]
# ## Define the nonlinear variational problem

# %%
# Set up nonlinear problem
problem = NonlinearProblem(Res, w, bcs, a)

# the global newton solver and params
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.convergence_criterion = "incremental"
solver.rtol = 1e-8
solver.atol = 1e-8
solver.max_it = 50
solver.report = True
solver.error_on_nonconvergence = False

#  The Krylov solver parameters.
ksp = solver.krylov_solver
opts = PETSc.Options()
option_prefix = ksp.getOptionsPrefix()
opts[f"{option_prefix}ksp_type"] = "preonly" # "preonly" works equally well
opts[f"{option_prefix}pc_type"] = "lu" # do not use 'gamg' pre-conditioner
opts[f"{option_prefix}pc_factor_mat_solver_type"] = "mumps"
opts[f"{option_prefix}ksp_max_it"] = 30
ksp.setFromOptions()


# %% [markdown]
# ##  Start calculation loop

# %%
# Give the step a descriptive name
step = "Indent"

# Variables for storing time history
totSteps = 1000000
timeHist0 = np.zeros(shape=[totSteps])
timeHist1 = np.zeros(shape=[totSteps]) 
timeHist2 = np.zeros(shape=[totSteps]) 
timeHist3 = np.zeros(shape=[totSteps]) 

#Iinitialize a counter for reporting data
ii=0

#  Set up temporary "helper" functions and expressions 
#  for updating the internal variables.
#
# For the Cv tensors:
Cv_1_temp = Function(V3)
Cv_1_expr = Expression(Cv_1, V3.element.interpolation_points())
#
Cv_2_temp = Function(V3)
Cv_2_expr = Expression(Cv_2, V3.element.interpolation_points())
#
Cv_3_temp = Function(V3)
Cv_3_expr = Expression(Cv_3, V3.element.interpolation_points())
#
Cv_4_temp = Function(V3)
Cv_4_expr = Expression(Cv_4, V3.element.interpolation_points())
#
Cv_5_temp = Function(V3)
Cv_5_expr = Expression(Cv_5, V3.element.interpolation_points())
#
# and also for the velocity and acceleration.
v_temp = Function(V2)
a_temp = Function(V2)
#
v_expr = Expression(v_new,V2.element.interpolation_points())
a_expr = Expression(a_new,V2.element.interpolation_points())

# Write initial state to file
writeResults(t=0.0)    

# print a message for simulation startup
mprint("------------------------------------")
mprint("Simulation Start")
mprint("------------------------------------")
# Store start time 
startTime = datetime.now()

# Time-stepping solution procedure loop
while (round(t + dt, 9) <= Ttot):
     
    # # larger dt during the zero-pressure period
    # if t + dt > tRamp:
    #   dt = dt2
    #   dk.value = dt2
    # else:
    #   step = "Indent"
     
    # increment time
    t += dt 

    
    # update time variables in time-dependent BCs 
    x3_ind.value =  length + r_ind - dispRamp(t - float(alpha)*dt)
    
    # Solve the problem
    (iter, converged) = solver.solve(w)
        
    
    # Now we start the adaptive time-stepping and output storage procedure.
    #
    # First, we check if the newton solver actually converged.
    if converged: 

        # Collect results from MPI ghost processes
        w.x.scatter_forward()
        
        # Write output to file
        writeResults(t)

        # increment counter
        ii += 1
        
        # If the solver converged, we print the status of the solver, 
        # perform adaptive time-stepping updates, output results, and 
        # update degrees of freedom for the next step, w_old <- w.
        
        # print progress of calculation periodically
        if ii%1 == 0:      
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            mprint("Step: {} |   Increment: {} | Iterations: {}".format(step, ii, iter))
            mprint("Simulation Time: {} s | dt: {} s".format(round(t,2), round(dt, 3)))
            mprint()   
        
        # Iteration-based adaptive time-stepping
        #
        # If the newton solver takes 4 or fewer iterations, 
        # increase the time step by a factor of 1.5:
        if ((iter<=4) and (dt<25)):
            dt = 1.5*dt
            dk.value = dt
        # If the newton solver takes 6 or more iterations, 
        # decrease the time step by a factor of 2:
        elif iter>=6:
            dt = dt/2
            dk.value = dt
        # otherwise the newton solver took 4-5 iterations,
        # in which case leave the time step alone.

        # Store time history variables at this time  
        timeHist0[ii] = t # current time
        #
        timeHist1[ii] =  dispRamp(t) # time history of indentation depth applied
        #
        try:
            timeHist2[ii] = w.sub(0).sub(2).eval([0.0, 0.0, length],colliding_cells[0])[0] # time history of displacement
        except:
            pass
        #
        try:
            timeHist3[ii] = domain.comm.gather(fem.assemble_scalar(rxnForce))[0] # time history of reaction force
        except:
            pass
        
        # update internal variables 
        #
        # interpolate the values of the internal variables into their "temp" functions
        Cv_1_temp.interpolate(Cv_1_expr)
        Cv_2_temp.interpolate(Cv_2_expr)
        Cv_3_temp.interpolate(Cv_3_expr)
        Cv_4_temp.interpolate(Cv_4_expr)
        Cv_5_temp.interpolate(Cv_5_expr)
        #
        v_temp.interpolate(v_expr)
        a_temp.interpolate(a_expr)

        # Update DOFs for next step
        w_old.x.array[:] = w.x.array
        
        # update the old values of internal variables for next step
        Cv_1_old.x.array[:] = Cv_1_temp.x.array[:]
        Cv_2_old.x.array[:] = Cv_2_temp.x.array[:]
        Cv_3_old.x.array[:] = Cv_3_temp.x.array[:]
        Cv_4_old.x.array[:] = Cv_4_temp.x.array[:]
        Cv_5_old.x.array[:] = Cv_5_temp.x.array[:]
        #
        v_old.x.array[:] = v_temp.x.array[:]
        a_old.x.array[:] = a_temp.x.array[:]
    
    # If solver doesn't converge we have to back up in time,
    # cut the size of the time step, and try solving again.
    else: # not(converged)
        
        # first, we back up in time
        # ( to un-do the current time step )
        t = t - float(dk)
        
        # Then, we cut back on the time step we're attempting.
        # (by a factor of 2)
        dt = dt/2
        dk.value = dt 

        # Re-set the DOFs to their value before the failed step.
        w.x.array[:] = w_old.x.array

        # Set up nonlinear problem
        problem = NonlinearProblem(Res, w, bcs, a)

        # The global newton solver and params
        solver = NewtonSolver(MPI.COMM_WORLD, problem)
        solver.convergence_criterion = "incremental"
        solver.rtol = 1e-8
        solver.atol = 1e-8
        solver.max_it = 50
        solver.report = True
        solver.error_on_nonconvergence = False
     
# close the output file.
file_results.close()
         
# End analysis
mprint("-----------------------------------------")
mprint("End computation")                 
# Report elapsed real time for the analysis
endTime = datetime.now()
elapseTime = endTime - startTime
mprint("------------------------------------------")
mprint("Elapsed real time:  {}".format(elapseTime))
mprint("------------------------------------------")



