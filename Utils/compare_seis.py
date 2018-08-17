#!/usr/bin/env python

import sys
import struct
import math
import os

#tolerance = 0.013
tolerance = 0.0011

#header in file 1
filename1 = sys.argv[1]
filename2 = sys.argv[2]

fp1 = open(filename1, 'r')
fp2 = open(filename2, 'r')

header_size=56
fp1.seek(header_size, os.SEEK_SET)
fp2.seek(header_size, os.SEEK_SET)

index = 0

num_comps = 2
nt = 2000

for c in range(0, num_comps):
	for t in range(0, nt):
		data1 = fp1.read(4)
		data2 = fp2.read(4)
		float1 = float(struct.unpack('f', data1)[0])
		float2 = float(struct.unpack('f', data2)[0])
		if float1!=0.0:
			if math.fabs(float1)<1.0:
				if math.fabs(float1-float2)>tolerance:
					print "Index %d: %f and %f differ by more than %f." % (index, float1, float2, tolerance)
		                        sys.exit(1)
			else:
				if math.fabs(float1-float2)/float1>tolerance:
					print "Index %d: %f and %f differ by more than %f%%." % (index, float1, float2, tolerance*100)
					sys.exit(1)
		index += 1
print "Seismograms agree."
