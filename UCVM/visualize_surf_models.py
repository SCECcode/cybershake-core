#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os

if len(sys.argv)<5:
	print "Usage: %s <surface file> <model coords file> <nx> <ny>" % (sys.argv[0])
	sys.exit(1)

surf_file = sys.argv[1]
model_coords_file = sys.argv[2]
nx = int(sys.argv[3])
ny = int(sys.argv[4])

#get list of points
coords = []
print "Reading model coods file."
with open(model_coords_file, "r") as fp_in:
	for line in fp_in:
		pieces = line.split()
		coords.append([float(pieces[0]), float(pieces[1])])

cvms_pts = []
usgs_pts = []
oned_pts = []
cca_pts = []
cvmsi_pts = []
other_pts = []

print "Reading surface file."
with open(surf_file, "r") as fp_in:
	for i in range(0, ny):
		if i%100==0:
			print "Line %d of %d." % (i, ny)
		line = fp_in.readline()
		pieces = line.split()
		for j in range(0, len(pieces)):
			if pieces[j]=="cvms5":
				cvms_pts.append(coords[i*nx+j])
			elif pieces[j]=="cencal":
				usgs_pts.append(coords[i*nx+j])
			elif pieces[j]=="1d":
				oned_pts.append(coords[i*nx+j])
			elif pieces[j]=="cca":
				cca_pts.append(coords[i*nx+j])
			elif pieces[j]=="cvmsi":
				cvmsi_pts.append(coords[i*nx+j])
			elif pieces[j]=="-1":
				other_pts.append(coords[i*nx+j])
			else:
				print "Label %s not recognized, aborting." % pieces[j]
				sys.exit(1)


clf()
plot([a[0] for a in oned_pts], [a[1] for a in oned_pts], '.', color="black", label="1d")
plot([a[0] for a in usgs_pts], [a[1] for a in usgs_pts], '.', color="green", label="usgs")
plot([a[0] for a in cvms_pts], [a[1] for a in cvms_pts], '.', color="red", label="cvms5")
plot([a[0] for a in cca_pts], [a[1] for a in cca_pts], '.', color="orange", label="cca")
plot([a[0] for a in cvmsi_pts], [a[1] for a in cvmsi_pts], '.', color="blue", label="cvmsi")
plot([a[0] for a in other_pts], [a[1] for a in other_pts], '.', color="gray", label="other")


legend()
xlim(-125, -113)
ylim(30, 40)
xlabel("Latitude")
ylabel("Longitude")
gcf().set_size_inches(8,8)
savefig("surface_plot.png", format="png")

