#!/usr/bin/env python3

'''This script will produce a profile plot at the simulation site.'''

import sys
import os

if len(sys.argv)<5:
	print("Usage: %s <media file> <gridout file> <fdloc file> <output profile filename>" % sys.argv[0])
	sys.exit(1)

media_file = sys.argv[1]
gridout_file = sys.argv[2]
fdloc_file = sys.argv[3]
output_filename = sys.argv[4]

#Get grid dimensions and spacing from gridout file
with open(gridout_file, 'r') as fp_in:
	data = fp_in.readlines()
	index = 0
	while index<len(data):
		line = data[index]
		if line[0:2]=="nx":
			xdim = int(line.strip().split("=")[1])
			spacing = float(data[index+1].split()[2])
			index += xdim
		elif line[0:2]=="ny":
			ydim = int(line.strip().split("=")[1])
			index += ydim
		elif line[0:2]=="nz":
			zdim = int(line.strip().split("=")[1])
			break
		else:
			index += 1
	fp_in.close()

with open(fdloc_file, 'r') as fp_in:
	data = fp_in.readlines()
	(x_coord, y_coord) = data[0].split()
	fp_in.close()

cmd = "%s/plot_profile.py %s %d %d %d %f %s %s %s" % (os.path.dirname(sys.argv[0]), media_file, xdim, ydim, zdim, spacing, x_coord, y_coord, output_filename)
sys.exit(os.system(cmd))
