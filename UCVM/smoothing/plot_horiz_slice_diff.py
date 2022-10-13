#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os
import struct
import math
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits import basemap
from mpl_toolkits.basemap import cm
from operator import itemgetter

'''Plots percent difference between two horizontal slices.'''

if len(sys.argv)<9:
	print("Usage: %s <velocity file 1> <velocity file 2> <model coords file> <nx> <ny> <nz> <z-depth> <decimation> <output file> [vp|vs|rho]" % sys.argv[0])
	sys.exit(1)


velocity_file_1 = sys.argv[1]
velocity_file_2 = sys.argv[2]
model_coords_file = sys.argv[3]
nx = int(sys.argv[4])
ny = int(sys.argv[5])
nz = int(sys.argv[6])
zslice = int(sys.argv[7])
decimation = int(sys.argv[8])
output_file = sys.argv[9]

comp = "vs"
if len(sys.argv)==11:
	comp = sys.argv[10]

x_dim = int(math.ceil(float(nx)/float(decimation)))
y_dim = int(math.ceil(float(ny)/float(decimation)))

coords = []
print("Reading model coods file.")
with open(model_coords_file, "r") as fp_in:
	data = fp_in.readlines()
	for y in range(0, ny, decimation):
		print("%d of %d y-points" % (y, ny))
		for x in range(0, nx, decimation):
			line_index = x + y*nx
			line = data[line_index]
			pieces = line.split()
			coords.append([float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])])

slice_floats = []

print("Reading velocity files.")
fp_in_1 = open(velocity_file_1, "rb")
fp_in_2 = open(velocity_file_2, "rb")
#Seek to z-slice
offset = zslice*12*nx*ny
fp_in_1.seek(offset)
fp_in_2.seek(offset)
slice_data_1 = fp_in_1.read(nx*ny*3*4)
slice_data_2 = fp_in_2.read(nx*ny*3*4)
#Plot % Vs diff
diffs = []
for i in range(0, nx, decimation):
    print("%d of %d x-slices" % (i, nx))
    offset = i*ny*3*4
    slice_floats_1 = (struct.unpack("%df" % 3*ny, slice_data_1[offset:offset+ny*3*4]))
    slice_floats_2 = (struct.unpack("%df" % 3*ny, slice_data_2[offset:offset+ny*3*4]))
    if comp=="vp":
        for f in range(0, len(slice_floats_1), 3*decimation):
            diffs.append(100.0*(slice_floats_2[f]-slice_floats_1[f])/slice_floats_1[f])
    elif comp=="vs":
        for f in range(1, len(slice_floats_1), 3*decimation):
            diffs.append(100.0*(slice_floats_2[f]-slice_floats_1[f])/slice_floats_1[f])
    elif comp=="rho":
        for f in range(2, len(slice_floats_1), 3*decimation):
            diffs.append(100.0*(slice_floats_2[f]-slice_floats_1[f])/slice_floats_1[f])
	
ordered_coords = sorted(coords, key=itemgetter(2,3))

#Plotting code taken from UCVM horizontal slice
#BOUNDS = [0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]
#BOUNDS = [0, 3000.0]
BOUNDS = [-100.0, -50.0, -20.0, -10.0, -5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
#TICKS = [0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]
#TICKS = [0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0]
TICKS = [-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]

#colormap = cm.RdBu
#Use new color map from SRL paper
colormap = basemap.cm.GMT_seis
norm = mcolors.Normalize(vmin=BOUNDS[0],vmax=BOUNDS[len(BOUNDS) - 1])


#statewide: 30, 42, -130, -113
#so cal: 31, 36, -121, -115

min_lat=31
max_lat=36
min_lon=-121
max_lon=-115

m = basemap.Basemap(projection='cyl', llcrnrlat=min_lat, urcrnrlat=max_lat, llcrnrlon=min_lon, urcrnrlon=max_lon, resolution='f', anchor='C')
#m = basemap.Basemap(projection='cyl', llcrnrlat=34, urcrnrlat=36, llcrnrlon=-123, urcrnrlon=-121, resolution='f', anchor='C')

lat_ticks = np.arange(min_lat, max_lat, 2)
#lat_ticks = np.arange(34, 36, 1)
lon_ticks = np.arange(min_lon, max_lon, 5)
#lon_ticks = np.arange(-123, -121, 1)

m.drawparallels(lat_ticks, linewidth=1.0, labels=[1,0,0,0])
m.drawmeridians(lon_ticks, linewidth=1.0, labels=[0,0,0,1])
m.drawstates()
m.drawcountries()

#pcolormesh needs 2d arrays
x_coords = np.array([a[0] for a in ordered_coords])
#x_2d = x_coords.reshape((nx, ny))
x_2d = x_coords.reshape((x_dim, y_dim))
y_coords = np.array([a[1] for a in ordered_coords])
#y_2d = y_coords.reshape((nx, ny))
y_2d = y_coords.reshape((x_dim, y_dim))
#vs_2d = np.array(vs).reshape((nx, ny))
vel_2d = np.array(diffs).reshape((x_dim, y_dim))

#m.pcolormesh(x_2d, y_2d, vs_2d, cmap=colormap, norm=norm, shading='flat', edgecolors='face')
#plt.pcolormesh([a[0] for a in ordered_coords], [a[1] for a in ordered_coords], vs, cmap=colormap, norm=norm)
#t = m.transform_scalar(vs_2d, x_coords, y_coords, len(x_coords), len(y_coords))
#img = m.imshow(vs_2d, cmap=colormap, norm=norm, interpolation='nearest')
m.scatter(x_coords, y_coords, c=diffs, cmap=colormap, norm=norm, s=1, edgecolor='', marker='o')

#Adding lines to show where the cross-sections are
'''
x_cross_section_lats = []
x_cross_section_lons = []
y_cross_section_lats = []
y_cross_section_lons = []
x_vals = [1000, 2000]
y_vals = [1740, 3480, 5220]
for i in range(0, len(x_vals)):
	x_cross_section_lats.append([])
	x_cross_section_lons.append([])
for i in range(0, len(y_vals)):
	y_cross_section_lats.append([])
	y_cross_section_lons.append([])

for c in ordered_coords:
	for i in range(0, len(x_vals)):
		if c[2]==x_vals[i]:
			x_cross_section_lons[i].append(c[0])
			x_cross_section_lats[i].append(c[1])
	for j in range(0, len(y_vals)):
		if c[3]==y_vals[j]:
			y_cross_section_lons[j].append(c[0])
                        y_cross_section_lats[j].append(c[1])

for i in range(0, len(x_vals)):
	print "%d: (%f, %f) to (%f, %f)" % (x_vals[i], x_cross_section_lats[i][0], x_cross_section_lons[i][0], x_cross_section_lats[i][-1], x_cross_section_lons[i][-1])
	m.plot(x_cross_section_lons[i], x_cross_section_lats[i], latlon=True, label="X=%d" % x_vals[i])

for i in range(0, len(y_vals)):
        print "%d: (%f, %f) to (%f, %f)" % (y_vals[i], y_cross_section_lats[i][0], y_cross_section_lons[i][0], y_cross_section_lats[i][-1], y_cross_section_lons[i][-1])
	m.plot(y_cross_section_lons[i], y_cross_section_lats[i], latlon=True, label="Y=%d" % y_vals[i])
'''
m.drawcoastlines()

cax = plt.axes([0.125, 0.05, 0.775, 0.02])
cbar = plt.colorbar(cax=cax, orientation='horizontal', ticks=TICKS)
gcf().set_size_inches(10, 8)
plt.savefig(output_file, type="png")
