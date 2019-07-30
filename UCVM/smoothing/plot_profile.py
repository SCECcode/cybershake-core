#!/usr/bin/env python

import sys
import os
import struct
import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *

if len(sys.argv)<8:
	print "Usage: %s <velocity file 1> <nx> <ny> <nz> <x-coord 1> <y-coord 1> <output file>" % sys.argv[0]
	sys.exit(1)

velocity_file = sys.argv[1]
nx = int(sys.argv[2])
ny = int(sys.argv[3])
nz = int(sys.argv[4])
x_coord = int(sys.argv[5])
y_coord = int(sys.argv[6])
output_filename = sys.argv[7]

print "Reading velocity file."

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
	depths.append(-0.1*i)


clf()
labels = ["Vp (km/s)", "Vs (km/s)", "rho (kg/m3)"]
for i in range(0, 3):
	plot(profile[i], depths, label=labels[i])
ylabel("Depth (km)")
xlim(0, 8000)
ylim(-0.1*nz, 0)
legend(loc="lower left")
#grid(b=True, which='major')
#minorticks_on()
#grid(b=True, which='minor')
gcf().set_size_inches(8, 10)
savefig(output_filename, type="png")
#Second plot of top 2 km
ylim(-2.0, 0)
xlim(0, 6000)
legend(loc='upper right')
savefig("zoomed_%s" % output_filename, type="png")



