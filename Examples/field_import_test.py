from PyPATools.field import Field
import os
import numpy as np
from numpy import meshgrid
import matplotlib.pyplot as plt

mydebug = True

# Test OPERA 3D field import
folder = r"/mnt/d/MIT Dropbox/Daniel Winklehner/Projects/IsoDAR/Cyclotron BFields"

bfield = Field.from_file(os.path.join(folder, "isodar2012_plugField_compl.table"),
                         label="Test Cyclotron B-Field 3D",
                         debug=mydebug,
                         scaling=-1.0)

print(bfield)
print(bfield(np.array([0.0, 0.0, 0.0])))

x = np.linspace(-0.14, 0.14, 100)
y = np.linspace(-0.14, 0.14, 100)


mesh_x, mesh_y = meshgrid(x, y, indexing='ij')
points = np.column_stack([mesh_x.flatten(), mesh_y.flatten(), np.zeros(10000)])

bz = bfield(points)[:, 2]

plt.contour(x, y, bz.reshape([100, 100]), 40)
plt.colorbar()
plt.show()

x = np.zeros(200)
y = np.zeros(200)
z = np.linspace(-1.5, 1.5, 200)

points = np.column_stack([x, y, z])

bz = bfield(points)[:, 2]

plt.plot(z, bz)
plt.show()
