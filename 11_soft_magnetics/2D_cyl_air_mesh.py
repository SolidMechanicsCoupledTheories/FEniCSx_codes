"""
2D Dogbone mesh

"""
"""
Step 1. Create a mesh with pygmsh 
"""
import pygmsh

# Resolution parameters
resolution = 0.5 # coarse element size
ratio = 20       # how much finer the fine elements at the corners are
air_res = 5          # fair-field air element size

# Geometry parameters
ld = 0.80
w = 5.9   # lateral width
l = ld*2*w # axial length
r = 0.2 # radius of corner fillet
s = 50  # side length of air domain

# Create an empty geometry, and initialize it
geometry = pygmsh.geo.Geometry()
# Create  a model to add data to
model = geometry.__enter__()

# Add points 
points = [model.add_point((0, 0), mesh_size=resolution),         # 0
          model.add_point((0, l/2), mesh_size=resolution),       # 1
          model.add_point((w-r, l/2), mesh_size=resolution/ratio),   # 2
          model.add_point((w, l/2-r), mesh_size=resolution/ratio),   # 3
          model.add_point((w, 0), mesh_size=resolution),         # 4
          model.add_point((w, -l/2+r), mesh_size=resolution/ratio),  # 5
          model.add_point((w-r, -l/2), mesh_size=resolution/ratio),  # 6
          model.add_point((0,-l/2), mesh_size=resolution),       # 7
          model.add_point((w-r, l/2-r), mesh_size=resolution/ratio), # 8
          model.add_point((w-r, -l/2+r), mesh_size=resolution/ratio), # 9
          model.add_point((0,s), mesh_size=air_res),      # 10
          model.add_point((s,s), mesh_size=air_res),      # 11
          model.add_point((s,-s), mesh_size=air_res),      # 12
          model.add_point((0,-s), mesh_size=air_res)]      # 13

# Add circles and lines 
line1 = model.add_line(points[0], points[1])
line2 = model.add_line(points[1], points[2])
circle1 = model.add_circle_arc(points[2], points[8], points[3])
line3 = model.add_line(points[3], points[4])
line4 = model.add_line(points[4], points[5])
circle2 = model.add_circle_arc(points[5], points[9], points[6])
line5 = model.add_line(points[6], points[7])
line6 = model.add_line(points[7], points[0])
#
line7 = model.add_line(points[10], points[1])
line8 = model.add_line(points[11], points[10])
line9 = model.add_line(points[12], points[11])
line10 = model.add_line(points[13], points[12])
line11 = model.add_line(points[7], points[13])

# Create a line_loop and plane_surface for meshing
lines_loop = model.add_curve_loop([line1, line2, circle1, line3, line4, circle2, line5, line6])
plane_surface = model.add_plane_surface(lines_loop)

air_loop = model.add_curve_loop([line2, circle1, line3, line4, circle2, line5,\
                                 line11, line10, line9, line8, line7])
plane_surface2 = model.add_plane_surface(air_loop)

# Call gmsh before adding physical entities
model.synchronize()

# The final step before mesh generation is to mark the domain  
# and the different boundaries. Give these entities names so that 
# they can be identified in gmsh 
model.add_physical(line8, "Top")
model.add_physical(line10, "Bottom")
model.add_physical([plane_surface], "MRE")
model.add_physical([plane_surface2], "Air")
model.add_physical([line2, circle1, line3, line4, circle2, line5 ], "Boundary")
model.add_physical([line7, line1, line6, line11], "Left")
model.add_physical(line9, "Right")

"""
In Fenics and Paraview these geometrical entities are nubered as:
   Top      = 1
   Bottom   = 2
   MRE      = 3
   Air      = 4
   Boundary = 5
   Left     = 6
   Right    = 7
   
"""

"""
Step 2.  Write the mesh to file using gmsh
"""
import gmsh
geometry.generate_mesh(dim=2)
gmsh.write("meshes/cyl_air.msh")
gmsh.clear()
geometry.__exit__()

# Now that we have saved the mesh to a `msh` file, we would like
# to convert it to a format that interfaces with dolfin/fenics. 

"""
Step 3. convert the mesh to.xdmf  format using meshio
"""
import meshio
mesh_from_file = meshio.read("meshes/cyl_air.msh")


"""
Step. 4  Extract cells and boundary data.

Now that we have created the mesh, we need to extract the cells 
and physical data. We need to create a separate file for the 
facets (lines),  which we will use when we define boundary 
conditions in  Fenics. We do this  with the following convenience 
function. Note that as we would like a  2 dimensional mesh, we need to 
remove the z-values in the mesh coordinates, if any.
"""
import numpy
def create_mesh(mesh, cell_type, prune_z=False):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    points = mesh.points[:,:2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells},\
                           cell_data={"name_to_read":[cell_data]})
    return out_mesh

"""
Step 5.
With this function in hand, we can save the facet line mesh 
and the domain triangle  mesh in `XDMF` format 
"""

line_mesh = create_mesh(mesh_from_file, "line", prune_z=True)
meshio.write("meshes/facet_cyl_air_0p80.xdmf", line_mesh)

triangle_mesh = create_mesh(mesh_from_file, "triangle", prune_z=True)
meshio.write("meshes/cyl_air_0p80.xdmf", triangle_mesh)


