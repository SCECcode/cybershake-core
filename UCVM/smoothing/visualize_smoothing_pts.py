#!/usr/bin/env python3

import matplotlib
matplotlib.use("AGG")
from pylab import *
import sys
import os

if len(sys.argv)<2:
	print("Usage: %s <smoothing points>" % (sys.argv[0]))
	sys.exit(1)

smoothing_file = sys.argv[1]

lon = []
lat = []
with open(smoothing_file, "r") as fp_in:
	for line in fp_in:
		pieces = line.split()
		lon.append(float(pieces[2]))
		lat.append(float(pieces[3]))

clf()
plot(lon, lat, '.')
xlim(-130, -110)
ylim(30, 42)
xlabel("Latitude")
ylabel("Longitude")
gcf().set_size_inches(8,8)
savefig("smoothing_points.png", format="png")

