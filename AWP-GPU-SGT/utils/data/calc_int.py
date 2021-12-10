#!/usr/bin/env python

import sys
import os

if len(sys.argv)<3:
	print("Usage: %s <source file> <dt>" % sys.argv[0])
	sys.exit(1)

dt = float(sys.argv[2])
with open(sys.argv[1], "r") as fp_in:
	data = fp_in.readlines()
	fp_in.close()
	sum = 0.0
	for line in data:
		pieces = line.split(",")
		sum += float(pieces[0])/1.0e15*dt

print(sum)
