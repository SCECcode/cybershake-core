#!/usr/bin/env python

import sys
import os
from operator import itemgetter


if len(sys.argv)<3:
	print "Usage: %s <rwg cordfile> <awp cordfile>" % sys.argv[0]
	sys.exit(1)

fp_in = open(sys.argv[1], "r")
fp_out = open(sys.argv[2], "w")

data = fp_in.readlines()
fp_in.close()

points = set()

for i in range(5, len(data), 1):
	line = data[i]
	pieces = line.split()
	x = int(pieces[0])
	y = int(pieces[1])
	z = int(pieces[2])
	points.add("%d %d %d\n" % (y+1, x+1, min([z+1, 200])))

p_list = []
for p in points:
	pieces = p.split()
	p_list.append([int(pieces[0]), int(pieces[1]), int(pieces[2])])
p_list = sorted(p_list, key=itemgetter(1,0,2))
fp_out.write("%d\n" % len(p_list))
for line in p_list:
	for p in line:
		fp_out.write("%d " % p)
	fp_out.write("\n")
fp_out.flush()
fp_out.close()
