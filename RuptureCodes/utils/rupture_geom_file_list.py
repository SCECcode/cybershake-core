#!/usr/bin/env python

import sys
import os
import re

if len(sys.argv)<7:
	print("Usage: %s <erf_id> <rup_var_scenario_id> <prefix> <pool> <rup geom file dir> <output file>" % sys.argv[0])
	print("Example: %s 35 4 gsiftp://hpc-login2.usc.edu hpc /home/rcf-104/CyberShake2007/ruptures/RuptureVariations_35_V3_2 hpc-e35-rv4.rls" % sys.argv[0])
	sys.exit(1)

#lfn, pfn, pool

erf_id = int(sys.argv[1])
rv_id = int(sys.argv[2])
prefix = sys.argv[3]
pool = sys.argv[4]
top_dir = sys.argv[5]
fp_out = open(sys.argv[6], "w")

regex = "\d+_\d+.txt$"
p = re.compile(regex)
index = 0
tot = len(os.listdir(top_dir))
for source_dir in os.listdir(top_dir):
	if os.path.isdir(os.path.join(top_dir, source_dir)):
		index += 1
		print("Source directory %s." % source_dir)
		if index%10==0:
			print("Source directory %d of %d." % (index, tot))
		for rup_dir in os.listdir(os.path.join(top_dir, source_dir)):
			if not os.path.isdir(os.path.join(top_dir, source_dir, rup_dir)):
				continue
			for entry in os.listdir(os.path.join(top_dir, source_dir, rup_dir)):
				if p.match(entry):
					print(entry)
					fp_out.write("e%d_rv%d_%s %s%s pool=%s\n" % (erf_id, rv_id, entry, prefix, os.path.join(top_dir, source_dir, rup_dir, entry), pool))
					break

fp_out.flush()
fp_out.close()

