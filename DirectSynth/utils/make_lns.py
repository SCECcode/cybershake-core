#!/usr/bin/env python

import sys
import os

rupture_file = sys.argv[1]
rupture_root = sys.argv[2]

with open(rupture_file, "r") as fp_in:
	data = fp_in.readlines()
	fp_in.close()
	for line in data[1:]:
		lfn = line.split()[0]
		#e36_rv6_39_5.txt
		pieces = lfn.split("_")
		src = int(pieces[2])
		rup = int(pieces[3].split(".")[0])
		target = "%s/%d/%d/%d_%d.txt" % (rupture_root, src, rup, src, rup)
		#print(target)
		os.symlink(target, lfn)
