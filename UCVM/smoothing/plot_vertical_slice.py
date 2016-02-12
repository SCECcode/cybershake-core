#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os
import struct
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits import basemap
from operator import itemgetter

if len(sys.argv)<9:
	print "Usage: %s <velocity file> <model coords file> <nx> <ny> <nz> <x or y> <row or column value> <output file>" % sys.argv[0]
	sys.exit(1)


velocity_file = sys.argv[1]
model_coords_file = sys.argv[2]
nx = int(sys.argv[3])
ny = int(sys.argv[4])
nz = int(sys.argv[5])
comp = sys.argv[6]
column = int(sys.argv[7])
output_file = sys.argv[8]

if comp!='x' and comp!='y':
	print "Component must be x or y."
	sys.exit(1)

if comp=='x':
	not_comp = 'y'
elif comp=='y':
	not_comp = 'x'

coords = []
print "Reading model coods file."
with open(model_coords_file, "r") as fp_in:
        for line in fp_in:
                pieces = line.split()
		if comp=='x' and int(pieces[3])==column:
	                coords.append([float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])])
		elif comp=='y' and int(pieces[2])==column:
			coords.append([float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])])

vs = []
x_values = []
y_values = []

GRID_SPACING = 0.4

print "Reading velocity file."
with open(velocity_file, "rb") as fp_in:
	for z in range(0, nz):
		for c in coords:
			offset = 4*3*(z*nx*ny + c[2]*ny + c[3])
			fp_in.seek(offset)
			data_str = fp_in.read(12)
			data = struct.unpack("3f", data_str)
			vs.append(data[1])
			if comp=='x':
				x_values.append(c[2]*GRID_SPACING)
			elif comp=='y':
				x_values.append(c[3]*GRID_SPACING)
			y_values.append(-1.0*z*GRID_SPACING)

#Plotting code taken from UCVM horizontal slice
BOUNDS = [0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]
TICKS = [0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]

colormap = cm.RdBu
norm = mcolors.Normalize(vmin=BOUNDS[0],vmax=BOUNDS[len(BOUNDS) - 1])

#plt.pcolormesh([a[0] for a in ordered_coords], [a[1] for a in ordered_coords], vs, cmap=colormap, norm=norm)
#t = m.transform_scalar(vs_2d, x_coords, y_coords, len(x_coords), len(y_coords))
#img = m.imshow(vs_2d, cmap=colormap, norm=norm, interpolation='nearest')
plt.figure(figsize=(12,5.5))
plt.clf()
plt.scatter(x_values, y_values, c=vs, cmap=colormap, norm=norm, marker='s', s=2, edgecolors='')
plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.15)
plt.xlabel("Distance (km) along %s=%d slice" % (not_comp, column))
plt.ylabel("Depth (km)")
plt.ylim(-80, 0)
if comp=='x':
	plt.xlim(0, nx*GRID_SPACING)
elif comp=='y':
	plt.xlim(0, ny*GRID_SPACING)
cax = plt.axes([0.125, 0.05, 0.775, 0.02])
cbar = plt.colorbar(cax=cax, orientation='horizontal', ticks=TICKS)

plt.savefig(output_file, type="png")
