#!/usr/bin/env python

'''This script combines the 3 steps required to make a smoothed mesh, which are:

1) Run determine_surface_model to create a file describing which model was used at each point along the surface.  From this, we can 
find the model interfaces.
2) Run determine_smoothing_points.py to translate the model surface file into a series of grid coordinates which require smoothing.
3) Run smooth (or smooth_mpi), which takes the list of points, and the mesh, and applies the smoothing.
'''

import sys
import os
import optparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

CS_PATH = config.getProperty("CS_PATH")

def run_det_surf_model(gridout, coords, models):
	'''Usage: %s <gridout file> <model coords file> <comma-separated list of velocity models> <output file>'''
	output_file = "surface_models.txt"
	cmd = "%s/UCVM/bin/determine_surface_model %s %s %s %s" % (CS_PATH, gridout, coords, models, output_file)
	print cmd
	exitcode = os.system(cmd)
	if (exitcode!=0):
		print "Error generating surface model, aborting."
		sys.exit(2)
	return output_file

def run_det_smooth_pts(surf_file, coords, nx, ny, smooth_dist):
	'''Usage: %s  <surf file> <model coords file> <nx> <ny> <smoothing dist> <output file>'''
	output_file = "smoothing_pts.txt"
	cmd = "%s/UCVM/smoothing/determine_smoothing_points.py %s %s %d %d %s %s" % (CS_PATH, surf_file, coords, nx, ny, smooth_dist, output_file)
	print cmd
	exitcode = os.system(cmd)
	if (exitcode!=0):
		print "Error generating list of smoothing points, aborting."
		sys.exit(3)
	return output_file

def run_smooth(input_mesh, smooth_pts, nx, ny, nz, smooth_dist, output_mesh):
	'''Usage: %s <AWP-format mesh in> <list of smoothing pts> <RWG nx> <RWG ny> <nz> <smoothing range in pts> <velocity mesh out>'''
	MPI_CMD = config.getProperty("MPI_CMD")
	#16 PPN
	np = int(os.environ["PBS_NUM_NODES"])*16
        MPI_CMD = "%s -n %d -N 16 -S 4" % (MPI_CMD, np)
	cmd = "%s %s/UCVM/smoothing/smooth_mpi %s %s %d %d %d %s %s" % (MPI_CMD, CS_PATH, input_mesh, smooth_pts, nx, ny, nz, smooth_dist, output_mesh)
	print cmd
	exitcode = os.system(cmd)
	if (exitcode!=0):
		print "Error performing smoothing, aborting."
		sys.exit(4)


parser = optparse.OptionParser()
parser.add_option("--gridout", dest="gridout", action="store", help="gridout file")
parser.add_option("--coords", dest="coords", action="store", help="coords file")
parser.add_option("--models", dest="modelString", action="store", help="comma-separated list of velocity models")
parser.add_option("--smoothing-dist", dest="smoothing_dist", action="store", type="int", help="Number of grid points to smooth over.  About 10km of grid points is a good starting place.")
parser.add_option("--mesh", dest="mesh", action="store", help="AWP-format velocity mesh to smooth")
parser.add_option("--mesh-out", dest="mesh_out", action="store", help="Output smoothed mesh")

(option, args) = parser.parse_args()

gridout = option.gridout
coords = option.coords
models = option.modelString
smoothing_dist = option.smoothing_dist
mesh_in = option.mesh
mesh_out = option.mesh_out

if (gridout==None) or (coords==None) or (mesh_in==None) or (mesh_out==None):
	print "Files gridout, model_coords, input velocity mesh, and output velocity mesh must be specified, aborting."
	sys.exit(1)

if (smoothing_dist==None):
	print "Smoothing distance in mesh points must be provided, aborting."
	sys.exit(1)

if (models==None):
	print "Comma-separated list of velocity models must be provided, aborting."
	sys.exit(1)

#Determine nx, ny, nz from gridout file
with open(gridout, "r") as fp_in:
	data = fp_in.readlines()
	fp_in.close()
	nx = int((data[1].split("="))[1])
	ny = int((data[1+nx+2].split("="))[1])
	nz = int((data[1+nx+2+ny+2].split("="))[1])

#Set LD_LIBRARY_PATH to pick up UCVM libraries
os.environ["LD_LIBRARY_PATH"] = "/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/lib/euclid3/lib:/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/lib/proj-4/lib:/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/model/cvms426/lib:/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/model/cencal/lib:/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/model/cvms5/lib:/projects/sciteam/bahm/CyberShake/software/UCVM/ucvm-15.10.0/model/cca/lib:%s" % (os.environ["LD_LIBRARY_PATH"])

surf_model_file = run_det_surf_model(gridout, coords, models)
smooth_pts_file = run_det_smooth_pts(surf_model_file, coords, nx, ny, smoothing_dist)
run_smooth(mesh_in, smooth_pts_file, nx, ny, nz, smoothing_dist, mesh_out)

