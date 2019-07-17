import fenics as fn
import numpy as np
from matplotlib import pyplot as plt
from spyral_inflector import *

x = TwoDeeField()

#
# gmsh_str = """
# SetFactory("OpenCASCADE");
# Mesh.CharacteristicLengthMax = 0.0025;
# dee_len = 0.2;
# h_gap = 0.01;
# v_gap = 0.05;
# dee_thk = 0.02;
#
# // Top dee electrode
# //Point(1) = {-(h_gap / 2 + dee_len), v_gap / 2 + dee_thk, 0.0};
# Point(2) = {-(h_gap / 2 + dee_len), v_gap / 2, 0.0};
# Point(3) = {-(h_gap / 2), v_gap / 2, 0.0};
# Point(4) = {-(h_gap / 2), v_gap / 2 + dee_thk, 0.0};
#
# //Line(1) = {1, 2};
# Line(2) = {2, 3};
# Line(3) = {3, 4};
# //Line(4) = {4, 1};
# //Line Loop(1) = {1:4};
#
# // Bottom dee electrode
# //Point(5) = {-(h_gap / 2 + dee_len), -(v_gap / 2 + dee_thk), 0.0};
# Point(6) = {-(h_gap / 2 + dee_len), -(v_gap / 2), 0.0};
# Point(7) = {-(h_gap / 2), -(v_gap / 2), 0.0};
# Point(8) = {-(h_gap / 2), -(v_gap / 2 + dee_thk), 0.0};
#
# //Line(5) = {5, 6};
# Line(6) = {6, 7};
# Line(7) = {7, 8};
# //Line(8) = {8, 5};
# //Line Loop(2) = {5:8};
#
# // Top dummy dee electrode
# //Point(9) = {h_gap / 2 + dee_len, v_gap / 2 + dee_thk, 0.0};
# Point(10) = {h_gap / 2 + dee_len, v_gap / 2, 0.0};
# Point(11) = {h_gap / 2, v_gap / 2, 0.0};
# Point(12) = {h_gap / 2, v_gap / 2 + dee_thk, 0.0};
#
# //Line(9) = {9, 10};
# Line(10) = {10, 11};
# Line(11) = {11, 12};
# //Line(12) = {12, 9};
# //Line Loop(3) = {9:12};
#
# // Bottom dummy dee electrode
# //Point(13) = {h_gap / 2 + dee_len, -(v_gap / 2 + dee_thk), 0.0};
# Point(14) = {h_gap / 2 + dee_len, -(v_gap / 2), 0.0};
# Point(15) = {h_gap / 2, -(v_gap / 2), 0.0};
# Point(16) = {h_gap / 2, -(v_gap / 2 + dee_thk), 0.0};
#
# //Line(13) = {13, 14};
# Line(14) = {14, 15};
# Line(15) = {15, 16};
# //Line(16) = {16, 13};
# //Line Loop(4) = {13:16};
#
# // Neumann boundary lines, clockwise
#
# Line(100) = {6, 2}; // left
# Line(101) = {16, 8}; // bottom
# Line(102) = {4, 12}; // top
# Line(103) = {10, 14}; // right
#
# Line Loop(1) = {2, 3, 102, -11, -10, 103, 14, 15, 101, -7, -6, 100};
# Surface(1) = {1};
#
# Physical Curve(1) = {2, 3, 6, 7}; // Dee
# Physical Curve(2) = {10, 11, 14, 15}; // Dummy Dee
# Physical Curve(3) = {102, 101}; // Top/Bottom boundary
# Physical Curve(4) = {100, 103}; // Left/Right boundary
#
# Physical Surface(1) = {1}; // Domain
#
# """
# TEMP_DIR = "/home/philip/work"
#
# with open('/home/philip/work/2dee_struct.geo', 'w') as f:
#     f.write(gmsh_str)
#
# import os
# import meshio
#
# os.system('gmsh -2 '+TEMP_DIR+'/2dee_struct.geo -format msh2 -v 0 -o '+TEMP_DIR+'/master_geometry.msh')
# # os.system('dolfin-convert master_geometry.msh master_geometry.xml')
# msh = meshio.read(TEMP_DIR + "/master_geometry.msh")
# # print(msh)
# # print(msh.cells)
# # meshio.write(TEMP_DIR + "/master_geometry.xdmf", msh)
#
# meshio.write_points_cells(TEMP_DIR + "/master_geometry_markers.xdmf",
#                           msh.points,
#                           {"triangle": msh.cells["triangle"]},
#                           cell_data={"triangle": {"gmsh:physical": msh.cell_data["triangle"]["gmsh:physical"]}})
#
# meshio.write_points_cells(TEMP_DIR + "/master_geometry_boundaries.xdmf",
#                           msh.points,
#                           {"line": msh.cells["line"]},
#                           cell_data={"line": {"gmsh:physical": msh.cell_data["line"]["gmsh:physical"]}})
#
# mesh = fn.Mesh()
# fn.XDMFFile(TEMP_DIR + "/master_geometry_markers.xdmf").read(mesh)
#
# markers = fn.MeshFunction("size_t", mesh, mesh.topology().dim())
# fn.XDMFFile(TEMP_DIR + "/master_geometry_markers.xdmf").read(markers, "gmsh:physical")
#
# _boundaries = fn.MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
# fn.XDMFFile(TEMP_DIR + "/master_geometry_boundaries.xdmf").read(_boundaries, "gmsh:physical")
# boundaries = fn.MeshFunction("size_t", mesh, _boundaries)
#
# # full_mesh = [mesh, markers, boundaries]
#
# dx = fn.Measure('dx', domain=mesh, subdomain_data=markers)
# V = fn.FunctionSpace(mesh, 'P', 1)
#
# # Build boundary conditions
# # TODO: Get the #'s from the meshing stuff -PW
#
# bcs = [fn.DirichletBC(V, fn.Constant(1.0), boundaries, 1),
#        fn.DirichletBC(V, fn.Constant(0), boundaries, 2)]
#
#
# # Test and trial functions
# u = fn.TrialFunction(V)
# v = fn.TestFunction(V)
#
# a = fn.dot(fn.grad(u), fn.grad(v)) * dx
# L = fn.Constant('0.0') * v * dx
#
# u = fn.Function(V)
#
# fn.solve(a == L, u, bcs, solver_parameters={"linear_solver": "cg", "preconditioner": "ilu"})
#
# print("Done!", flush=True)
#
# potentialFile = fn.File(TEMP_DIR + '/potential.pvd')
# potentialFile << u
#
# meshfile = fn.File(TEMP_DIR + '/mesh.pvd')
# meshfile << mesh
#
# fenics_field = fn.project(-fn.grad(u), solver_type='cg', preconditioner_type='ilu')
# electric_field = FenicsField(fenics_field)
#
# xr = np.linspace(-0.1, 0.1, 100)
# ef = []
# for x in xr:
#     ef.append(electric_field(np.array([x, 0.0, 0.0]))[0])
#
# plt.plot(xr, ef)
# plt.plot(xr, np.max(ef)*np.exp(-650*xr**2))
# plt.show()