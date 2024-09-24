#!/usr/bin/env python3

import math
import os
import sys
from pyproj import Proj

if len(sys.argv)<10:
	print("Usage: %s <lat to search> <lon to search> <center lat> <center lon> <center X index> <center Y index> <rotation angle (degrees CCW)> <grid spacing in m> <utm zone>" % (sys.argv[0]))
	sys.exit(1)

point_lat = float(sys.argv[1])
point_lon = float(sys.argv[2])
center_lat = float(sys.argv[3])
center_lon = float(sys.argv[4])
center_x = int(sys.argv[5])
center_y = int(sys.argv[6])
rot_angle = float(sys.argv[7])
spacing = int(sys.argv[8])
zone = int(sys.argv[9])

rot_rad = rot_angle*math.pi/180
#construct inverse of rotational matrix
rot = [[math.cos(rot_rad), math.sin(rot_rad)], [-1.0*math.sin(rot_rad), math.cos(rot_rad)]]

#Determine UTM coords of point and center
proj = Proj(proj='utm', zone=zone, ellps='WGS84')
(point_east, point_north) = proj(point_lon, point_lat)
(center_east, center_north) = proj(center_lon, center_lat)

#Difference
diff_east = point_east - center_east
diff_north = point_north - center_north

#Distance along X/Y axes
x_dist = rot[0][0]*diff_east + rot[0][1]*diff_north
y_dist = rot[1][0]*diff_east + rot[1][1]*diff_north

#Grid offset along axes
x_grid_diff = x_dist/spacing
y_grid_diff = y_dist/spacing

#Closest grid point
point_x = center_x + x_grid_diff
point_y = center_y + y_grid_diff

print "%f %f %f %f" % (point_lon, point_lat, point_x, point_y)

