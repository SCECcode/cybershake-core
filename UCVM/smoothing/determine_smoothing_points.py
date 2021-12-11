#!/usr/bin/env python3

import sys
import os
from operator import itemgetter

BEFORE = 0
ABOVE = 1

smoothing_points = set()

def add_points(y, x, direction):
        for i in range(x-1, x-SMOOTHING_DIST-1, -1):
                if i<0:
                        continue
                else:
                        if coords[y][i] not in smoothing_points:
                            smoothing_points.add(coords[y][i])
        if direction==BEFORE:
		#Do 1 more earlier
                i = x-SMOOTHING_DIST-1
                if i>=0:
                    smoothing_points.add(coords[y][i])			
        for i in range(x, x+SMOOTHING_DIST+1):
                if i>=nx: 
                        continue
                else:   
                        if coords[y][i] not in smoothing_points:
                                smoothing_points.add(coords[y][i])
        for j in range(y-1, y-SMOOTHING_DIST-1, -1):
                if j<0:
                        continue
                else:
                        if coords[j][x] not in smoothing_points:
                            smoothing_points.add(coords[j][x])
        if direction==ABOVE:
                j = y-SMOOTHING_DIST-1	
                if j>=0:
                    smoothing_points.add(coords[j][x])
        for j in range(y, y+SMOOTHING_DIST+1):
                if j>=ny:
                        continue
                else:
                        if coords[j][x] not in smoothing_points:
                                smoothing_points.add(coords[j][x])


#number of mesh points away from boundary to smooth
#175m spacing in the mesh
#10150m on either size
SMOOTHING_DIST = 58

if len(sys.argv)<7:
	print("Usage: %s <surf file> <model coords file> <nx> <ny> <smoothing dist> <output file>" % sys.argv[0])
	sys.exit(1)

surf_file = sys.argv[1]
model_coords = sys.argv[2]
nx = int(sys.argv[3])
ny = int(sys.argv[4])
SMOOTHING_DIST = int(sys.argv[5])
output_file = sys.argv[6]

#get list of points
coords = []
for i in range(0, ny):
	coords.append([])

print("Reading model coords file.")
with open(model_coords, "r") as fp_in:
        data = fp_in.readlines()
        for i in range(0, ny):
                for j in range(0, nx):
                        line = data[j+i*nx]
                        pieces = line.split()
                        coords[i].append((float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])))
        fp_in.close()

print("Reading surface file.")
with open(surf_file, "r") as fp_in:
	data = fp_in.readlines()
	line_before = data[0]
	pieces_before = line_before.split()
	for i in range(1, ny):
		line = data[i]
		pieces = line.split()
		#Check for boundaries
		for j in range(1, nx):
			#Check entry before and entry above
			#Look at just usgs/cvms boundary
			if pieces[j]!=pieces[j-1]:
				add_points(i, j, BEFORE)
			elif pieces[j]!=pieces_before[j]:
				add_points(i, j, ABOVE)
			'''
			if pieces[j]=="cencal" and pieces[j-1]=="1d":
				add_points(i, j, BEFORE)
			elif pieces[j]=="cencal" and pieces_before[j]=="1d":
				add_points(i, j, ABOVE)
			elif pieces[j]=="1d" and pieces[j-1]=="cencal":
				add_points(i, j, BEFORE)
                        elif pieces[j]=="1d" and pieces_before[j]=="cencal":
                                add_points(i, j, ABOVE)
			'''
		pieces_before = pieces
	fp_in.close()

print(len(smoothing_points))

#Write points to file
with open(output_file, "w") as fp_out:
	ordered_list = sorted(smoothing_points, key=itemgetter(2,3))
	for p in ordered_list:
		fp_out.write("%d %d %f %f\n" % (int(p[2]), int(p[3]), float(p[0]), float(p[1])))
	fp_out.flush()


