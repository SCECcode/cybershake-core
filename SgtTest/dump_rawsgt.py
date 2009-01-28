#!/usr/bin/env python

import sys

sys.path.append('..')

import config

if len(sys.argv)<3:
	print "Usage: %s <prefix> <nproc>\n" % sys.argv[0]
	sys.exit(1)

pref = sys.argv[1]
nproc = sys.argv[2]
cs_path = config.getProperty('CS_PATH')


cmd = './dump_rawsgt.csh %s %s %s' % (pref, nproc, cs_path)
exitcode = os.system(cmd)
sys.exit((exitcode >> 8) & 0xFF)

