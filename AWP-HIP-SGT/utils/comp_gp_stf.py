#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Feb 24 2020

@author: Zhifeng Hu <zhh076@ucsd.edu>
Input:  *.srf : srf source file from CyberShake
        surface grid : mesh grids, from UCVM

Output: gp_src.bin, GP-type source for IFAULT=5
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
from scipy import spatial
import os
import sys
from pathlib2 import Path
import pickle
import struct

stf_type = 3.0  # GP source

# Parameters. Users should define themselves
mx, my = 720, 720  # dimensions (x, y); 
nt_ref = 6000  # desired time steps of source 
dt = 0.005  # desired time interval
dz = 0.100  # grid spacing in km
theta_rot = 0.0  # rotation angle, positive for clockwise
block_top = 0  # Top the the mesh block, not needed for uniform mesh
# If the fault size is small, specifyting a rough box can speed up calculation
xoff_0, xoff_1 = 0, mx  # left, right boundary
yoff_0, yoff_1 = 0, my  # top, bottom boundary
kx, ky = xoff_1 - xoff_0, yoff_1 - yoff_0

# Read input file names
if len(sys.argv) < 3:
    print('Input parameters required: srf, surface grid!')
    sys.exit(-1)
srf_file = sys.argv[1]
surf_grid = sys.argv[2]

# Read surface grid locations, to search for subfault indices
try:
    kdtree = pickle.loads(Path('kdtree.pickle')).read_bytes()
except:
    grids = np.fromfile(surf_grid, dtype='float64').reshape(my, mx, 3)
    grids = np.reshape(grids[yoff_0 : yoff_1, xoff_0 : xoff_1, :2], (-1, 2))
    kdtree = spatial.cKDTree(grids)
    with open('kdtree.pickle', 'wb') as fid:
        pickle.dump(kdtree, fid, protocol=pickle.HIGHEST_PROTOCOL)
    del grids

f = open(srf_file,'r')
f.readline()
f.readline()
token = f.readline()
f.readline()
print(f.readline())
nx = int(token.split()[2])
nz = int(token.split()[3])

# Output files
res = np.zeros((nz, nx, 12), dtype='float32')
idx = np.zeros((nz, nx, 3), dtype='int32')

for j in range(nz):
    sys.stdout.write("\rreading subfault %d of %d" % (j+1, nz))
    sys.stdout.flush()
    last_pos = f.tell()
    for i in range(nx):
        nl1 = f.readline().split()
        nl2 = f.readline().split()
        iy, ix = np.unravel_index(kdtree.query([float(nl1[0]), float(nl1[1])])[1], (ky, kx))
        iz = int(round((float(nl1[2]) - block_top) / dz))
        stk = float(nl1[3]) - theta_rot
        dip = float(nl1[4])
        area = float(nl1[5]) / 1.e4  # cm^2 --> m^2
        tinit= float(nl1[6])
        dt = float(nl1[7])
        rake = float(nl2[0])
        slip = float(nl2[1]) / 1.e2  # cm --> m
        nt1 = int(nl2[2]) 
        nskip1 = int(np.ceil(nt1 / 6.))
        for l in range(nskip1):
            tmp = f.readline()
        trise = nt1 * dt
        idx[j, i, :] = np.array([ix + xoff_0 + 1, iy + yoff_0 + 1, iz + 1]) # 1-index
        res[j, i, :] = [stf_type, slip, tinit, trise, area, 0.0, stk, dip, rake, 0.0, 0.0, 0.0]
f.close()

# Reshape output arrays
idx = idx.reshape((-1, 3))
res = res.reshape((-1, 12))

# Write gp_src.bin
fmt = "<3i12f"
with open('gp_src.bin', 'wb') as f_out:
    for i in range(len(res)):
        f_out.write(struct.pack(fmt, *idx[i, :], *res[i, :])) 

# In `param.sh`, set INSRC gp_src.bin; IFAULT 5
try:
    os.symlink('gp_src.bin', 'gp_src.bin_0')
except FileExistsError:
    os.remove('gp_src.bin_0')
    os.symlink('gp_src.bin', 'gp_src.bin_0')
    pass
