#!/usr/bin/env python

import sys
import os
import struct
import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *

if len(sys.argv)<10:
	print "Usage: %s <velocity file 1> <model coords file> <nx> <ny> <nz> <x-coord 1> <y-coord 1> <output file>" % sys.argv[0]
	sys.exit(1)

velocity_file_1 = sys.argv[1]
model_coords_file = sys.argv[2]
nx = int(sys.argv[3])
ny = int(sys.argv[4])
nz = int(sys.argv[5])
x_coord_1 = int(sys.argv[6])
y_coord_1 = int(sys.argv[7])
x_coord_2 = int(sys.argv[8])
y_coord_2 = int(sys.argv[9])
output_filename = sys.argv[10]

'''
with open(model_coords_file, "r") as fp_in:
        for line in fp_in:
                pieces = line.split()
		if int(pieces[2])==x_coord and int(pieces[3])==y_coord:
			target_lon = float(pieces[0])
			target_lat = float(pieces[1])
'''
print "Reading velocity file 1."

profile_1 = []
profile_2 = []
for i in range(0, 3):
	profile_1.append([])
	profile_2.append([])

with open(velocity_file_1, "rb") as fp_in:
	for z in range(0, nz):
		offset = 4*3*(z*nx*ny + x_coord_1*ny + y_coord_1)
		fp_in.seek(offset, os.SEEK_SET)
		data_str = fp_in.read(12)
		data = struct.unpack("3f", data_str)
		for i in range(0, 3):
			profile_1[i].append(data[i])
		offset = 4*3*(z*nx*ny + x_coord_2*ny + y_coord_2)
                fp_in.seek(offset, os.SEEK_SET)
                data_str = fp_in.read(12)
                data = struct.unpack("3f", data_str)
                for i in range(0, 3):
                        profile_2[i].append(data[i])

'''
print "Reading velocity file 2."
profile_2 = []
for i in range(0, 3):
        profile_2.append([])

with open(velocity_file_2, "rb") as fp_in:
        for z in range(0, nz):
                offset = 4*3*(z*nx*ny + x_coord*ny + y_coord)
                fp_in.seek(offset, os.SEEK_SET)
                data_str = fp_in.read(12)
                data = struct.unpack("3f", data_str)
                for i in range(0, 3):
                        profile_2[i].append(data[i])
'''

depths = []
for i in range(0, nz):
	depths.append(-0.4*i)


clf()
subplots_adjust(left=0.05)
subplots_adjust(right=0.95)
subplots_adjust(top=0.92)
subplots_adjust(bottom=0.08)
subplot(131, title="Vp")
plot(profile_1[0], depths, label="x=%d, y=%d" % (x_coord_1, y_coord_1))
plot(profile_2[0], depths, label="x=%d, y=%d" % (x_coord_2, y_coord_2))
xlabel("Vp (km/s)")
ylabel("Depth (km)")
xlim(0, 8000)
legend(loc="lower left")
subplot(132, title="Vs")
plot(profile_1[1], depths, label="x=%d, y=%d" % (x_coord_1, y_coord_1))
plot(profile_2[1], depths, label="x=%d, y=%d" % (x_coord_2, y_coord_2))
xlabel("Vs (km/s)")
xlim(0, 5000)
#xlim(3000, 5000)
#ylim(-35, -25)
ylabel("Depth (km)")
#grid(b=True, which='major')
#minorticks_on()
#grid(b=True, which='minor')
legend(loc="lower left")
subplot(133, title="Density")
plot(profile_1[2], depths, label="x=%d, y=%d" % (x_coord_1, y_coord_1))
plot(profile_2[2], depths, label="x=%d, y=%d" % (x_coord_2, y_coord_2))
xlabel("Density")
ylabel("Depth (km)")
xlim(0, 5000)
legend(loc='lower left')
gcf().set_size_inches(13, 6)
savefig(output_filename, type="png")
