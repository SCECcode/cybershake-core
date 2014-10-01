#!/usr/bin/env python

import sys
import os
import config

def build_cordfile(site, rwg_cordfile, awp_cordfile, max_depth):
	cs_path = config.getProperty("CS_PATH")
	
	cmd = "%s/SgtHead/gen_awp_cordfile.py %s %s" % (cs_path, rwg_cordfile, awp_cordfile, max_depth)

	exitcode = os.system(cmd)
	return ((exitcode >> 8) & 0xFF)

