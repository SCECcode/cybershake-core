#!/usr/bin/env python

import sys
import os
import struct

if len(sys.argv)<3:
	print("Usage: %s <model coords> <surf grid file out>" % (sys.argv[0]))
	sys.exit(1)

model_coords = sys.argv[1]
surf_grid_file = sys.argv[2]

with open(model_coords, "r") as fp_in:
	with open(surf_grid_file, "wb") as fp_out:
		#surf grid file should be fast x, slow y
		#Since model_coords is fast y, slow x, can convert coordinates from RWG to AWP and use in same order
		data = fp_in.readlines()
		fp_in.close()
		for line in data:
			pieces = line.split()
			out_str = struct.pack("ddd", float(pieces[0]), float(pieces[1]), 0.0)
			fp_out.write(out_str)
		fp_out.flush()
		fp_out.close()


