#!/usr/bin/env python

import sys
import os
import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *

filename = sys.argv[1]

fp_in = open(filename, "r")
data = fp_in.readlines()
fp_in.close()
sgt1 = []
sgt2 = []

for i in range(0, 6):
	sgt1.append([])
	sgt2.append([])
	pieces = data[2*i].split()
	for piece in pieces:
		sgt1[i].append(float(piece))
        pieces = data[2*i+1].split()
        for piece in pieces:
                sgt2[i].append(float(piece))

	
min_y = 0.0
max_y = 0.0
for i in range(0, 6):
	if min(sgt1[i])<min_y:
                min_y = min(sgt1[i])
        if max(sgt1[i])>max_y:
                max_y = max(sgt1[i])
        if min(sgt2[i])<min_y:
                min_y = min(sgt2[i])
        if max(sgt2[i])>max_y:
                max_y = max(sgt2[i])


min_y = 1.1*min_y
max_y = 1.1*max_y

sgt1_ts = []
sgt2_ts = []

for j in range(0, 4000):
	sgt2_ts.append(0.05*j)
	if j%2==0:
		sgt1_ts.append(0.05*j)

clf()
components = ["XX", "YY", "ZZ", "XY", "XZ", "YZ"]
subplots_adjust(hspace=0.25, left=0.09, right=0.98, top=0.96, bottom=0.04)
for i in range(0, 3):
	for j in range(0, 2):
		index = 321 + i*2 + j
		subplot(index, title="%s" % components[i*2 + j])
		plot(sgt1_ts, sgt1[2*i+j], label="0.5Hz")
		plot(sgt2_ts, sgt2[2*i+j], label="1.0Hz")
		xlim(0, 200)
		ylim(min_y, max_y)
		legend()
gcf().set_size_inches(14, 7)
outfile = "%s.png" % (filename.split(".")[0])
savefig(outfile, format="png")

