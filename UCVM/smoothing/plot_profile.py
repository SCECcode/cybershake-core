#!/usr/bin/env python

import sys
import os
import struct
import matplotlib
matplotlib.use("AGG")
from pylab import *

if len(sys.argv)<9:
	print("Usage: %s <velocity file 1> <nx> <ny> <nz> <grid spacing in km> <x-coord 1> <y-coord 1> <output file>" % sys.argv[0])
	sys.exit(1)

velocity_file = sys.argv[1]
nx = int(sys.argv[2])
ny = int(sys.argv[3])
nz = int(sys.argv[4])
grid_spacing = float(sys.argv[5])
x_coord = int(sys.argv[6])
y_coord = int(sys.argv[7])
output_filename = sys.argv[8]
output_data_filename = "%s.data" % ".".join(output_filename.split(".")[:-1])

print("Reading velocity file.")

profile = []
for i in range(0, 3):
	profile.append([])

with open(velocity_file, "rb") as fp_in:
	for z in range(0, nz):
		offset = 4*3*(z*nx*ny + x_coord*ny + y_coord)
		fp_in.seek(offset, os.SEEK_SET)
		data_str = fp_in.read(12)
		data = struct.unpack("3f", data_str)
		for i in range(0, 3):
			profile[i].append(data[i])

depths = []
for i in range(0, nz):
	depths.append(-1.0*grid_spacing*i)

with open(output_data_filename, "w") as fp_out:
	fp_out.write("Depth(m),Vp(m/s),Vs(m/s),rho(kg/m3)\n")
	for i in range(0, nz):
		fp_out.write("%f,%f,%f,%f\n" % (depths[i], profile[0][i], profile[1][i], profile[2][i]))
	fp_out.flush()
	fp_out.close()


clf()
labels = ["Vp (m/s)", "Vs (m/s)", "rho (kg/m3)"]
for i in range(0, 3):
	plot(profile[i], depths, label=labels[i])
ylabel("Depth (km)")
xlim(0, 8000)
ymin = -3.0
#print(profile[0])
#ylim(-1.0*grid_spacing*nz, 0)
ylim(ymin, 0.0)
legend(loc="lower left")
#grid(b=True, which='major')
#minorticks_on()
#grid(b=True, which='minor')
gcf().set_size_inches(8, 10)
savefig(output_filename, format="png")
#Second plot of top 2 km
#ylim(-2.0, 0)
#xlim(0, 6000)
#legend(loc='upper right')
#savefig("zoomed_%s" % output_filename, format="png")



