#!/usr/bin/env python

import sys
import os

def build_cordfile(site, rwg_cordfile, awp_cordfile, max_depth):
	import config
	cs_path = config.getProperty("CS_PATH")
	
	cmd = "%s/SgtHead/gen_awp_cordfile.py %s %s %d" % (cs_path, rwg_cordfile, awp_cordfile, max_depth)
	print(cmd)
	exitcode = os.system(cmd)
	return ((exitcode >> 8) & 0xFF)

if __name__=="__main__":
	#Add config to path
	full_path = os.path.abspath(sys.argv[0])
	path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))
	sys.path.append(path_add)
	if len(sys.argv)<5:
		print("Usage: %s <site> <rwg cordfile> <awp cordfile> <max depth index, usually nz>" % (sys.argv[0]))
		sys.exit(1)
	build_cordfile(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
	
