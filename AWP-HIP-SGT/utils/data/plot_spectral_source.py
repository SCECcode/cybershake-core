#!/usr/bin/env python

import sys
import os
import numpy
import matplotlib
from pylab import *

if len(sys.argv)<2:
	print "Usage: %s <source file>" % sys.argv[0]

in_file = sys.argv[1]

fp_in = open(in_file, "r")
data = fp_in.readlines()
fp_in.close()

data_list = []

for line in data:
	pieces = line.split(",")
	data_list.append(float(pieces[0].strip()))


data_arr = numpy.array(data_list)
N = len(data_list)
dt = 0.01
fft = numpy.fft.fft(data_arr)
freqs = numpy.linspace(0.0, 2.0, num=2.0/dt)
max_ind = int(2.0/dt)
plot(freqs, abs(fft[0:max_ind]))
#matplotlib.pyplot.show()
gcf().set_size_inches(8,6)
outfile = "%s.png" % (in_file)
savefig(outfile, format="png")

