#!/usr/bin/env python

import sys
import os

if len(sys.argv)<3:
	print "Usage: %s <source file> <comp>" % sys.argv[0]
	sys.exit(1)

fp_in = open(sys.argv[1], "r")
fp_out = open("f%s_src" % sys.argv[2], "w")

data = fp_in.readlines()
fp_in.close()

pos = -1

if sys.argv[2]=="x":
	pos = 0
elif sys.argv[2]=="y":
	pos = 1
elif sys.argv[2]=="z":
	pos = 2

for i in range(2, len(data)):
	pieces = data[i].split()
	for p in pieces:
		for j in range(0, 3):
			if j==pos:
				fp_out.write("%e, " % (float(p)*1.0e15))
			else:
				fp_out.write("0.0, ")
		fp_out.write("0.0, 0.0, 0.0\n")
fp_out.flush()
fp_out.close()
