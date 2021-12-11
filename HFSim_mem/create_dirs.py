#!/usr/bin/env python

import sys
import os

if len(sys.argv)<2:
	print("Usage: %s <file with list of dirs>" % (sys.argv[0]))
	sys.exit(1)

for line in open(sys.argv[1], "r"):
	dirname = line.strip()
	if not os.path.exists(dirname):
		os.mkdir(line.strip())

