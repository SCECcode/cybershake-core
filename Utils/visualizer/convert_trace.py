#!/usr/bin/env python

'''Take CSV with <lat, lon> and convert into <X, Y> coordinates based on grid volume, for use with Paraview.'''

import sys
import os
import math

if len(sys.argv)<5:
	print "Usage: %s <lat, lon trace file> <model coords file> <grid spacing in km> <X, Y output csv>" % sys.argv[0]
	sys.exit(1)

grid_spacing = float(sys.argv[3])

coords_data = []
print "Reading coords data."
with open(sys.argv[2], "r") as fp_in:
	data = fp_in.readlines()
	for d in data:
		pieces = d.split()
		coords_data.append([float(pieces[0]), float(pieces[1]), int(pieces[2]), int(pieces[3])])
	fp_in.close()

PI = 3.1416
RAD = 6371

with open(sys.argv[1], "r") as fp_in:
	with open(sys.argv[4], "w") as fp_out:
		print "Reading trace data."
		data = fp_in.readlines()
		for m in range(0, len(data)):
			if (m%10)==0:
				print "%d of %d processed." % (m, len(data))
			d = data[m]
			nearest_index = [-1, -1, -1, -1]
			nearest_dists = [1000.0, 2000.0, 3000.0, 4000.0]
			pieces = d.split(',')
			#pieces[0] = lat, pieces[1] = lon
			lat = float(pieces[0])
			lon = float(pieces[1])
			#Find nearest points
			for k in range(0, len(coords_data)):
				c = coords_data[k]
				if math.fabs(lat-c[1])>.01:
					continue
				if math.fabs(lon-c[0])>.01:
					continue
				#calculate distance
				lat_dist = (lat-c[0])*PI/180.0*RAD
				med_lat = (lat+c[0])/2.0
				lon_dist = (lon-c[1])*PI/180.0*RAD*math.cos(med_lat*PI/180.0)
				dist = math.sqrt(lat_dist*lat_dist + lon_dist*lon_dist)
				if dist>2*grid_spacing:
					continue
				for i in range(0, 4):
					if dist<nearest_dists[i]:
						for j in range(2, -1, -1):
							nearest_dists[j+1] = nearest_dists[j]
							nearest_index[j+1] = nearest_index[j]
						nearest_dists[i] = dist
						nearest_indices[i] = k
						break
			#inverse-dist average the 4 to get the coords
			dist_tot = 0.0
			x_tot = 0.0
			y_tot = 0.0
			for i in range(0, len(nearest_index)):
				if nearest_index[i]==-1:
					break
				x_tot += 1.0/(nearest_dists[i])*coords_data[i][2]
				y_tot += 1.0/(nearest_dists[i])*coords_data[i][3]
				dist_tot += 1.0/(nearest_dists[i])
			if dist_tot==0.0 or x_tot==0.0 or y_tot==0.0:
				continue
			x = x_tot/dist_tot
			y = y_tot/dist_tot
			#Flip X and Y since model_coords is in RWG coordinates but we want AWP 
			fp_out.write("%.2f,%.2f\n" % (y, x))
		fp_out.flush()
		fp_out.close()
	fp_in.close()


