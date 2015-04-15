#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os
import struct

ref_seis = sys.argv[1]
test_seis = sys.argv[2]

ref_ts = []
test_ts = []
ref_data = []
test_data = []

for i in range(0, 2):
	ref_data.append([])
	test_data.append([])

ref_NT = 2000
test_NT = 2000
ref_DT = 0.1
test_DT = 0.1

for i in range(0, ref_NT):
	ref_ts.append(i*ref_DT)
for i in range(0, test_NT):
	test_ts.append(i*test_DT)

ref_in = open(ref_seis, "r")
#ignore header
ref_in.seek(56, os.SEEK_SET)
ref_str = ref_in.read(4*ref_NT)
ref_data[0] = struct.unpack("%df" % ref_NT, ref_str)
ref_str = ref_in.read(4*ref_NT)
ref_data[1] = struct.unpack("%df" % ref_NT, ref_str)
ref_in.close()

test_in = open(test_seis, "r")
test_in.seek(56, os.SEEK_SET)
test_str = test_in.read(4*test_NT)
test_data[0] = struct.unpack("%df" % test_NT, test_str)
test_str = test_in.read(4*test_NT)
test_data[1] = struct.unpack("%df" % test_NT, test_str)
test_in.close()

max_y = max([max(ref_data[0]), max(test_data[0]), max(ref_data[1]), max(test_data[1])])
max_y = 1.1*max_y

min_y = min([min(ref_data[0]), min(test_data[0]), min(ref_data[1]), min(test_data[1])])
min_y = 1.1*min_y

clf()
subplots_adjust(hspace=0.25, left=0.09, right=0.98, top=0.96, bottom=0.04)
subplot(211, title="X component")
plot(ref_ts, ref_data[0], label="1000m SGT, 1000m surf")
plot(test_ts, test_data[0], label="200m SGT, 200m surf")
xlim(0, 150)
ylim(min_y, max_y)
legend()
subplot(212, title="Y component")
plot(ref_ts, ref_data[1], label="1000m SGT, 1000m surf")
plot(test_ts, test_data[1], label="200m SGT, 200m surf")
xlim(0, 150)
ylim(min_y, max_y)
legend()
gcf().set_size_inches(14, 7)
outfile = "%s.png" % (ref_seis.split(".")[0])
savefig(outfile, format="png")

