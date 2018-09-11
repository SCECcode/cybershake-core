#!/usr/bin/env python2

import sys
import os

if len(sys.argv)<3:
	print "Usage: %s <checkpoint in> <checkpoint out>" % sys.argv[0]
	sys.exit(1)

file_in = sys.argv[1]
file_out = sys.argv[2]

entries = dict()

with open(file_in, "r") as fp_in:
	with open(file_out, "w") as fp_out:
		data = fp_in.readlines()
		for line in data:
			key = line.strip()
			if key not in entries:
				entries[key] = 1
		for key in entries:
			fp_out.write("%s\n", key)
		fp_out.flush()
		fp_out.close()
	fp_in.close()

