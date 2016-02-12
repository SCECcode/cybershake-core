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

if len(sys.argv)<8:
	print "Usage: %s <velocity file> <model coords file> <nx> <ny> <nz> <z-depth> <output file>" % sys.argv[0]
	sys.exit(1)


velocity_file = sys.argv[1]
model_coords_file = sys.argv[2]
nx = int(sys.argv[3])
ny = int(sys.argv[4])
nz = int(sys.argv[5])
zslice = int(sys.argv[6])
output_file = sys.argv[7]

coords = []
print "Reading model coods file."
with open(model_coords_file, "r") as fp_in:
        for line in fp_in:
                pieces = line.split()
                coords.append([float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])])

slice_floats = []
vs = []

print "Reading velocity file."
with open(velocity_file, "rb") as fp_in:
	#Seek to z-slice
	offset = zslice*12*nx*ny
	fp_in.seek(offset)
	slice_data = fp_in.read(nx*ny*3*4)
	#Plot Vs
	for i in range(0, nx):
		offset = i*ny*3*4
		slice_floats = (struct.unpack("%df" % 3*ny, slice_data[offset:offset+ny*3*4]))
		vs.extend(slice_floats[1::3])
	
ordered_coords = sorted(coords, key=itemgetter(2,3))

#Plotting code taken from UCVM horizontal slice
BOUNDS = [0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]
TICKS = [0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]

colormap = cm.RdBu
norm = mcolors.Normalize(vmin=BOUNDS[0],vmax=BOUNDS[len(BOUNDS) - 1])

m = basemap.Basemap(projection='cyl', llcrnrlat=30, urcrnrlat=42, llcrnrlon=-128, urcrnrlon=-112, resolution='f', anchor='C')
#m = basemap.Basemap(projection='cyl', llcrnrlat=34, urcrnrlat=36, llcrnrlon=-123, urcrnrlon=-121, resolution='f', anchor='C')

lat_ticks = np.arange(30, 42, 2)
#lat_ticks = np.arange(34, 36, 1)
lon_ticks = np.arange(-128, -112, 5)
#lon_ticks = np.arange(-123, -121, 1)

m.drawparallels(lat_ticks, linewidth=1.0, labels=[1,0,0,0])
m.drawmeridians(lon_ticks, linewidth=1.0, labels=[0,0,0,1])
m.drawstates()
m.drawcountries()

#pcolormesh needs 2d arrays
x_coords = np.array([a[0] for a in ordered_coords])
x_2d = x_coords.reshape((nx, ny))
y_coords = np.array([a[1] for a in ordered_coords])
y_2d = y_coords.reshape((nx, ny))
vs_2d = np.array(vs).reshape((nx, ny))

fp_out = open("surface_vel.txt", "w")
for i in range(0, len(vs)):
	fp_out.write("%f %f %f\n" %(y_coords[i], x_coords[i], vs[i]))
fp_out.flush()
fp_out.close()

#m.pcolormesh(x_2d, y_2d, vs_2d, cmap=colormap, norm=norm, shading='flat', edgecolors='face')
#plt.pcolormesh([a[0] for a in ordered_coords], [a[1] for a in ordered_coords], vs, cmap=colormap, norm=norm)
#t = m.transform_scalar(vs_2d, x_coords, y_coords, len(x_coords), len(y_coords))
#img = m.imshow(vs_2d, cmap=colormap, norm=norm, interpolation='nearest')
m.scatter(x_coords, y_coords, c=vs, cmap=colormap, norm=norm, s=1, edgecolor='', marker='.')

#Adding lines to show where the cross-sections are
#5 cross-sections, at X=183, 367, 550, 733, 917
'''
cross_section_lats = []
cross_section_lons = []
for i in range(0, 5):
	cross_section_lats.append([])
	cross_section_lons.append([])
x_vals = [183, 367, 550, 733, 917]
for c in ordered_coords:
	for i in range(0, len(x_vals)):
		if c[2]==x_vals[i]:
			cross_section_lons[i].append(c[0])
			cross_section_lats[i].append(c[1])

for i in range(0, len(x_vals)):
	m.plot(cross_section_lons[i], cross_section_lats[i], latlon=True, label="X=%d" % x_vals[i])
'''
m.drawcoastlines()

cax = plt.axes([0.125, 0.05, 0.775, 0.02])
cbar = plt.colorbar(cax=cax, orientation='horizontal', ticks=TICKS)

plt.savefig(output_file, type="png")
