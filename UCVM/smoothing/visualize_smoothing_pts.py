#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os

if len(sys.argv)<2:
	print "Usage: %s <smoothing points>" % (sys.argv[0])
	sys.exit(1)

smoothing_file = sys.argv[1]

x = []
y = []
with open(smoothing_file, "r") as fp_in:
	for line in fp_in:
		pieces = line.split()
		x.append(float(pieces[0]))
		y.append(float(pieces[1]))

clf()
plot(x, y, '.')
xlim(-130, -110)
ylim(30, 42)
xlabel("Latitude")
ylabel("Longitude")
gcf().set_size_inches(8,8)
savefig("smoothing_points.png", format="png")

