#!/usr/bin/env python

import sys
import os
from operator import itemgetter


if len(sys.argv)<4:
	print("Usage: %s <rwg cordfile> <awp cordfile> <max depth index>" % sys.argv[0])
	sys.exit(1)

fp_in = open(sys.argv[1], "r")
fp_out = open(sys.argv[2], "w")
max_depth_index = int(sys.argv[3])

data = fp_in.readlines()
fp_in.close()

points = set()

num_pts_in = int(data[4])

for i in range(5, len(data), 1):
	line = data[i]
	pieces = line.split()
	x = int(pieces[0])
	y = int(pieces[1])
	z = int(pieces[2])
        #Adjusted point_str to not add 1 to z-coordinate, since both AWP and RWG use z=1 to represent the free surface
	point_str = "%d %d %d\n" % (y+1, x+1, z)
	if point_str in points:
		print("Duplicate point entry %s" % point_str)
	points.add(point_str)

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

#Make sure both cordfiles have the same number of points
if num_pts_in!=len(p_list):
	print("Error: input file %s has %d points, but %d points were written to output file %s." % (sys.argv[1], num_pts_in, len(p_list), sys.argv[2]))
	sys.exit(2)

