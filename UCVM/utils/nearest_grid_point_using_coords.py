#!/usr/bin/env python3

import sys
import os
from pyproj import Proj
import math

if len(sys.argv)<5:
	print("Usage: %s <lat> <lon> <model coords file> <utm zone>" % (sys.argv[0]))
	sys.exit(1)

point_lat = float(sys.argv[1])
point_lon = float(sys.argv[2])
model_coords_file = sys.argv[3]
zone = int(sys.argv[4])

proj = Proj(proj='utm', zone=zone, ellps='WGS84')
(point_e, point_n) = proj(point_lon, point_lat)

closest_lat = -1.0
closest_lon = -1.0
closest_x = -1
closest_y = -1
closest_dist = 1.0e10

with open(model_coords_file, "r") as fp_in:
	data = fp_in.readlines()
	counter = 0
	for line in data:
		counter += 1
		if counter%100000==0:
			print("%d of %d" % (counter, len(data)))
		(test_lon_str, test_lat_str) = line.split()[0:2]
		test_lat = float(test_lat_str)
		if math.fabs(point_lat-test_lat)>0.005:
			continue
		test_lon = float(test_lon_str)
		if math.fabs(point_lon-test_lon)>0.005:
			continue
		(test_e, test_n) = proj(test_lon, test_lat)
		dist2 = ((test_e-point_e)*(test_e-point_e))+((test_n-point_n)*(test_n-point_n))
		if dist2<closest_dist:
			closest_dist = dist2
			#print test_lon, test_lat, test_e, point_e, test_n, point_n, closest_dist
			closest_lat = test_lat
			closest_lon = test_lon
			closest_x = int(line.split()[2])
			closest_y = int(line.split()[3])
	fp_in.close()

print("Closest point to (%f, %f) is (%f, %f), %f m away, with X=%d, Y=%d" % (point_lat, point_lon, closest_lat, closest_lon, math.sqrt(closest_dist), closest_x, closest_y))

