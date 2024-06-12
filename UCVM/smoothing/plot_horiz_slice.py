#!/usr/bin/env python2

import matplotlib
matplotlib.use("AGG")
from pylab import *
import sys
import os
import struct
import math
import numpy as np
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits import basemap
from mpl_toolkits.basemap import cm
from operator import itemgetter

if len(sys.argv)<8:
	print "Usage: %s <velocity file> <model coords file> <nx> <ny> <nz> <z-depth> <decimation> <output file> [vp|vs|rho]" % sys.argv[0]
	sys.exit(1)


velocity_file = sys.argv[1]
model_coords_file = sys.argv[2]
nx = int(sys.argv[3])
ny = int(sys.argv[4])
nz = int(sys.argv[5])
zslice = int(sys.argv[6])
decimation = int(sys.argv[7])
output_file = sys.argv[8]

comp = "vs"
if len(sys.argv)==10:
	comp = sys.argv[9]

x_dim = int(math.ceil(float(nx)/float(decimation)))
y_dim = int(math.ceil(float(ny)/float(decimation)))

coords = []
print "Reading model coods file."
with open(model_coords_file, "r") as fp_in:
	data = fp_in.readlines()
	for y in range(0, ny, decimation):
		print "%d of %d y-points" % (y, ny)
		for x in range(0, nx, decimation):
			line_index = x + y*nx
			line = data[line_index]
			pieces = line.split()
			coords.append([float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])])

slice_floats = []
vel_data = []

print "Reading velocity file."
with open(velocity_file, "rb") as fp_in:
	#Seek to z-slice
	offset = zslice*12*nx*ny
	fp_in.seek(offset)
	slice_data = fp_in.read(nx*ny*3*4)
	#Plot Vs
	for i in range(0, nx, decimation):
		print "%d of %d x-slices" % (i, nx)
		offset = i*ny*3*4
		slice_floats = (struct.unpack("%df" % 3*ny, slice_data[offset:offset+ny*3*4]))
		if comp=="vp":
	                vel_data.extend(slice_floats[0::3*decimation])
		elif comp=="vs":
	                vel_data.extend(slice_floats[1::3*decimation])
		elif comp=="rho":
			vel_data.extend(slice_floats[2::3*decimation])
	
ordered_coords = sorted(coords, key=itemgetter(2,3))

#Plotting code taken from UCVM horizontal slice
#BOUNDS = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]
#BOUNDS = [400.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0]
BOUNDS = [500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0]
#BOUNDS = [0, 3000.0]
#TICKS = [500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0]
#TICKS = [400, 500, 1000, 1500, 2000, 2500, 3000, 4000]
TICKS = [500, 1000, 1500, 2000, 2500, 3000, 4000, 5000]
#TICKS = [0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0, 3000.0]

#colormap = cm.RdBu
#Use new color map from SRL paper
#colormap = basemap.cm.GMT_seis
#Use color-blind-friendly color map
colormap = mcm.viridis
#norm = mcolors.Normalize(vmin=BOUNDS[0],vmax=BOUNDS[len(BOUNDS) - 1])
norm = mcolors.LogNorm(vmin=BOUNDS[0],vmax=BOUNDS[len(BOUNDS)-1])

#statewide: 30, 42, -130, -113
#so cal: 31, 36, -121, -115
#no cal: 33, 43, -127, -115

min_lat=31
max_lat=38
min_lon=-123
max_lon=-113

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
vel_2d = np.array(vel_data).reshape((x_dim, y_dim))

fp_out = open("surface_vel.txt", "w")
for i in range(0, len(vel_data)):
	fp_out.write("%f %f %f\n" %(y_coords[i], x_coords[i], vel_data[i]))
fp_out.flush()
fp_out.close()

#m.pcolormesh(x_2d, y_2d, vs_2d, cmap=colormap, norm=norm, shading='flat', edgecolors='face')
#plt.pcolormesh([a[0] for a in ordered_coords], [a[1] for a in ordered_coords], vs, cmap=colormap, norm=norm)
#t = m.transform_scalar(vs_2d, x_coords, y_coords, len(x_coords), len(y_coords))
#img = m.imshow(vs_2d, cmap=colormap, norm=norm, interpolation='nearest')
m.scatter(x_coords, y_coords, c=vel_data, cmap=colormap, norm=norm, s=1, edgecolor=None, marker='o')

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
cbar = plt.colorbar(cax=cax, orientation='horizontal', ticks=TICKS, format="%d")
#cbar.ax.set_xscale('log')
gcf().set_size_inches(10, 8)
plt.savefig(output_file, format="png")
