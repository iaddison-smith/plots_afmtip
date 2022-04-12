#!/usr/bin/env python
# coding: utf-8

import numpy as np
import meshio
import bempp.api

def generate_grid(face_file, vert_file):
    """
    Create bempp Grid object from .face and .vert files.
    """
    face = open(face_file, 'r').read()
    vert = open(vert_file, 'r').read()
    faces = np.vstack(np.char.split(face.split('\n')[0:-1]))[:,:3].astype(int) - 1
    verts = np.vstack(np.char.split(vert.split('\n')[0:-1]))[:,:3].astype(float)
    grid = bempp.api.Grid(verts.transpose(), faces.transpose())

    return grid

# Path: generate_vtk.py
face = 'zika_afm\\surf_gs1.0_noTER_split.face'  # .face file
vert = 'zika_afm\\surf_gs1.0_noTER_split.vert'  # .vert file
print('Loading Mesh')
grid = generate_grid(face, vert)   # use bempp to generate grid
print('Mesh loaded! \n ----------------------------------------------')

# read solution file
print('Loading Results')
result = np.loadtxt('zika_afm\\tip_charge-2.5_rad150\\zatsc386\\phi.txt') # 2 nm
result2 = np.loadtxt('zika_afm\\tip_charge-2.5_rad150\\zatsc1384\\phi.txt')
result = result - result2 # subtract the two results to obtain the delta potential
print('Results loaded! \n ----------------------------------------------')

# Extract potential and derivative
potential = result[:grid.number_of_elements]
derivative = result[grid.number_of_elements:2*grid.number_of_elements]

# Calculate boundary forces from potential and derivative
ep_ex = 80
ep_in = 4
kappa = 0.125
to_kcalmolA = 4*np.pi*332.0636817823836
to_pN = 69.467
f_db = to_pN*to_kcalmolA*-0.5*(ep_ex-ep_in)*(ep_in/ep_ex)*derivative**2*grid.volumes
f_ib = to_pN*to_kcalmolA*-0.5*kappa*ep_ex*potential**2*grid.volumes
f_bind = (np.sqrt((-2.1197685475579964)**2+(0.36092747252832624)**2+(0.00760169067022572)**2))

# generate points, cells, cell_data to create vtk using meshio
points = grid.vertices.T                      
cells = [("triangle", grid.elements.T.astype("int32"))]
cell_data = dict()
cell_data['potential'] = [potential]
cell_data['derivative'] = [derivative]
cell_data['f_db'] = [f_db/f_bind]
cell_data['f_ib'] = [f_ib/f_bind]

meshio.write_points_cells(
    'zika_vtk\\zika_condelta_porcentual.vtk',
    points=points,
    cells=cells,
    cell_data=cell_data
)

