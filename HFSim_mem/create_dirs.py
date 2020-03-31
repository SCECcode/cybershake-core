#!/usr/bin/env python3

import sys
import os

if len(sys.argv)<2:
	print("Usage: %s <file with list of dirs>" % (sys.argv[0]))
	sys.exit(1)

for line in open(sys.argv[1], "r"):
	os.mkdir(line.strip())


