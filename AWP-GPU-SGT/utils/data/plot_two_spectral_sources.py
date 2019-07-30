#!/usr/bin/env python

import sys
import os
import numpy
import matplotlib
from pylab import *

if len(sys.argv)<2:
	print "Usage: %s <source file 1> <source file 2>" % sys.argv[0]

in_file_1 = sys.argv[1]
in_file_2 = sys.argv[2]

fp_in_1 = open(in_file_1, "r")
data_1 = fp_in_1.readlines()
fp_in_1.close()

fp_in_2 = open(in_file_2, "r")
data_2 = fp_in_2.readlines()
fp_in_2.close()

data_list_1 = []
data_list_2 = []

for line in data_1:
	pieces = line.split(",")
	data_list_1.append(float(pieces[0].strip()))

for line in data_2:
        pieces = line.split(",")
        data_list_2.append(float(pieces[0].strip()))


data_arr_1 = numpy.array(data_list_1)
N = len(data_list_1)
dt = 0.005
MAX_FREQ = 3.0
fft_1 = numpy.fft.fft(data_arr_1)
freqs = numpy.linspace(0.0, MAX_FREQ, num=MAX_FREQ/dt)
max_ind = int(MAX_FREQ/dt)
data_arr_2 = numpy.array(data_list_2)
fft_2 = numpy.fft.fft(data_arr_2)
plot(freqs, abs(fft_1[0:max_ind]))
plot(freqs, abs(fft_2[0:max_ind]))
xlim(0, MAX_FREQ)
#matplotlib.pyplot.show()
gcf().set_size_inches(8,6)
outfile = "%s_comp_%s.png" % (in_file_1, in_file_2)
savefig(outfile, format="png")

