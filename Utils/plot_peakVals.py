#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os
import struct

ref_psa = sys.argv[1]
test_psa = sys.argv[2]
test_2_psa = sys.argv[3]

half_hz_periods = [10.0, 9.5, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.8, 4.6, 4.4, 4.2, 4.0, 3.8, 3.6, 3.4, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0]
one_hz_periods = [10.0, 9.5, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5, 6.0, 5.5, 5.0, 4.8, 4.6, 4.4, 4.2, 4.0, 3.8, 3.6, 3.4, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 10.0/6.0, 10.0/7.0, 10.0/8.0, 10.0/9.0, 1.0]
half_hz_freq = []
one_hz_freq = []
half_hz_data = []
half_hz_data_2 = []
one_hz_data = []

for p in half_hz_periods:
	half_hz_freq.append(1.0/p)

for p in one_hz_periods:
	one_hz_freq.append(1.0/p)

ref_in = open(ref_psa, "r")
#ignore header
ref_in.seek(56, os.SEEK_SET)
for i in range(0, len(half_hz_freq)):
	ref_str = ref_in.read(4)
	ref_data = struct.unpack("f", ref_str)
	half_hz_data.append(ref_data)	
ref_in.close()

test_in = open(test_psa, "r")
test_in.seek(56, os.SEEK_SET)
for i in range(0, len(half_hz_freq)):
	test_str = test_in.read(4)
	test_data = struct.unpack("f", test_str)
	half_hz_data_2.append(test_data)
test_in.close()

test_in_2 = open(test_2_psa, "r")
test_in_2.seek(56, os.SEEK_SET)
for i in range(0, len(one_hz_freq)):
        test_2_str = test_in_2.read(4)
        test_2_data = struct.unpack("f", test_2_str)
        one_hz_data.append(test_2_data)
test_in_2.close()


max_y = max([max(half_hz_data), max(half_hz_data_2), max(one_hz_data)])
#max_y = max([max(half_hz_data), max(half_hz_data_2)])
max_y = 1.1*max_y[0]
min_y = 0

clf()
plot(half_hz_freq, half_hz_data, label="0.5Hz, ERF35")
plot(half_hz_freq, half_hz_data_2, label="0.5Hz, ERF36")
plot(one_hz_freq, one_hz_data, label="1.0Hz, ERF36")
xlim(0.1, 1.0)
ylim(min_y, max_y)
xscale('log')
legend(loc='upper left')
gcf().set_size_inches(6, 6)
outfile = "%s.png" % (ref_psa.split(".")[0])
savefig(outfile, format="png")

