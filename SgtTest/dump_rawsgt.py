#!/usr/bin/env python

import sys
import os

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

if len(sys.argv)<3:
	print "Usage: %s <prefix> <nproc>\n" % sys.argv[0]
	sys.exit(1)

pref = sys.argv[1]
nproc = sys.argv[2]
cs_path = config.getProperty('CS_PATH')


cmd = '%s/dump_rawsgt.csh %s %s %s' % (sys.path[0], pref, nproc, cs_path)
exitcode = os.system(cmd)
sys.exit((exitcode >> 8) & 0xFF)

