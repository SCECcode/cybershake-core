#!/usr/bin/env python3

'''
This script reads in a list of surface points in a CyberShake model_coords file, and queries for the surface elevation in the UCVM DEM.
The 1d[SCEC] model is queried, since it's very fast.
'''

import sys
import os

from ucvm.src.framework.ucvm import UCVM
import ucvm.src.shared.constants
from ucvm.src.shared.properties import Point, SeismicData

def main(argv):
	if len(argv)<2:
		print("Usage: %s <model_coords file> <output_file>" % argv[0])
		sys.exit(1)

	coord_file = argv[1]
	output_file = argv[2]
	print("Parsing coords file.")
	coords = parse_coords(coord_file)
	print("Querying UCVM.")
	data = query_points(coords)
	print("Writing output.")
	write_output(data, output_file)


def parse_coords(coord_file):
	coords = []
	with open(coord_file, "r") as fp_in:
		data = fp_in.readlines()
		for line in data:
			pieces = line.split()
			coords.append([float(pieces[0]), float(pieces[1])])
	return coords

def query_points(coords_list):
	sd_objs = []
	for c in coords_list:
		sd_objs.append(SeismicData(Point(c[0], c[1], 0.0)))
	model_str = "1d[SCEC]"
	UCVM.query(sd_objs, model_str)
	return sd_objs
	

def write_output(data, output_file):
	with open(output_file, "w") as fp_out:
		for d in data:
			fp_out.write("%f %f %f\n" % (d.original_point.x_value, d.original_point.y_value, d.elevation_properties.elevation))
		fp_out.flush()
		fp_out.close()
	
if __name__=="__main__":
	main(sys.argv)
