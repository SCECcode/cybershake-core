#!/usr/bin/env python

import sys
import os
import numpy
import matplotlib
from pylab import *

if len(sys.argv)<4 or len(sys.argv)%2 != 0:
	print "Usage: %s <source file 1> <label 1> <source file 2> <label 2> ... <source file N> <label N> <output file>" % sys.argv[0]
	sys.exit(1)

filenames = []
labels = []
for f in range(1, len(sys.argv)-1, 2):
	filenames.append(sys.argv[f])
	labels.append(sys.argv[f+1])

outfile = sys.argv[-1]

#Assume all sources have the same length
MAXT = 200.0

for j in range(0, len(filenames)):
	f = filenames[j]
	fp_in = open(f, "r")
	data = fp_in.readlines()
	fp_in.close()
	data_vals = []
	for line in data:
		pieces = line.split(",")
	        data_vals.append(float(pieces[0].strip()))
	nt = len(data_vals)
	dt = MAXT/nt
	print "Using dt=%f, nt=%d" % (dt, nt)
	timesteps = [dt*i for i in range(0, nt)]
	plot(timesteps, data_vals, label=labels[j])

legend(loc="upper right")
xlim(0, MAXT/20.0)
gcf().set_size_inches(8,6)
savefig(outfile, format="png")

