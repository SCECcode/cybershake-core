#!/usr/bin/env python3

import sys
import os

if len(sys.argv)<3:
	print("Usage: %s <seis 0> <seis 1> ... <seis N> <output seis name>" % sys.argv[0])
	sys.exit(1)

seis_out = sys.argv[-1]
seis_in = sys.argv[1:-1]

cmd = "cat %s > %s" % (' '.join(seis_in), seis_out)
exitcode = os.system(cmd)
if exitcode!=0:
	print("Command '%s' failed." % (cmd))
sys.exit(exitcode)
