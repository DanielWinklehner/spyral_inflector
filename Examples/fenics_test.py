import fenics as fn
import os
import sys
import meshio


os.system('gmsh -3 working_example.geo -format msh2 -o working_example.msh')
# mesh = meshio.read('test_geo.msh')


os.system('dolfin-convert working_example.msh working_example.xml')
#
# os.system('gmsh -3 test_geo.geo -format msh2 -o test_geo.msh')
# os.system('dolfin-convert test_geo.msh test_geo.xml')
#

"""
Physical Volume(1) = {1};
anode_boundary[] = Boundary { Volume{ 1 }; };
N_anode = #anode_boundary[];
For i In {0:N_anode-1}
    Physical Surface (i+1) = { anode_boundary[i] };
EndFor

Cylinder(2) = {0, 0, -0.13, 0, 0, 0.16, 0.12, 2*Pi};
Physical Volume(2) = {2};
vacuum_boundary[] = Boundary { Volume{ 2 }; };
N_vacuum = #vacuum_boundary[];
For j In {0:N_vacuum-1}
    Physical Surface (N_anode +1 + j) = { vacuum_boundary[j] };
EndFor

Physical Volume(3) = {1001};
cathode_boundary[] = Boundary { Volume{ 1001 }; };
N_cathode = #cathode_boundary[];
For k In {0:N_cathode-1}
    Physical Surface (1000+k) = { cathode_boundary[k] };
EndFor

Delete{ Volume{1, 2, 1001}; }
Volume (1) = {2, 1, 3};
"""

mesh = fn.Mesh("working_example.xml")
markers = fn.MeshFunction('size_t', mesh, 'working_example_physical_region.xml')
boundaries = fn.MeshFunction('size_t', mesh, 'working_example_facet_region.xml')
dx = fn.Measure('dx', domain=mesh, subdomain_data=markers)
V = fn.FunctionSpace(mesh, 'P', 1)

bcs = []
for i in range(1, 258):
    bcs.append(fn.DirichletBC(V, fn.Constant(1E4), boundaries, i))
for i in range(258, 261):
    bcs.append(fn.DirichletBC(V, fn.Constant(0.0), boundaries, i))
for i in range(1000, 1257):
    bcs.append(fn.DirichletBC(V, fn.Constant(-1E4), boundaries, i))

u = fn.TrialFunction(V)
v = fn.TestFunction(V)
a = fn.dot(fn.grad(u), fn.grad(v)) * fn.dx
L = fn.Constant('0') * v * fn.dx
u = fn.Function(V)
fn.solve(a == L, u, bcs)

electric_field = -fn.grad(u)

potentialFile = fn.File('output/potential.pvd')
potentialFile << u

vtkfile = fn.File('output/e_field.pvd')
vtkfile << fn.project(electric_field)