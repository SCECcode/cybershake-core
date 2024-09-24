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
import math
from timeit import default_timer

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

CS_PATH = config.getProperty("CS_PATH")
MPI_CMD = config.getProperty("MPI_CMD")

def run_det_surf_model(gridout, coords, models):
	'''Usage: %s <gridout file> <model coords file> <comma-separated list of velocity models> <output file>'''
	output_file = "surface_models_mpi.txt"
	if MPI_CMD=="aprun":
		np = int(os.environ["PBS_NUM_NODES"])*8
		prefix = "%s -n %d -N 8 -S 4" % (MPI_CMD, np)
	elif MPI_CMD=="jsrun":
		np = int(os.environ["LSB_DJOB_NUMPROC"])-1
		#-a 1: 1 MPI tasks/core
		#-c 1: 1 core per resource set
		#-r 42: 42 resource sets per node
		prefix = "%s -a 1 -c 1 -r 42 -n %d" % (MPI_CMD, np)
		#Add export of PROJ_LIB
		prefix = "export PROJ_LIB=/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/lib/proj-4/share/proj; %s" % prefix
	elif MPI_CMD=="srun":
		num_nodes = int(os.environ["SLURM_JOB_NUM_NODES"])
		proc_per_node = int(os.environ["SLURM_CPUS_ON_NODE"])
		np = num_nodes * proc_per_node
		prefix = "%s -N %d -n %d -c 1" % (MPI_CMD, num_nodes, np)
	cmd = "/bin/date; %s %s/UCVM/bin/determine_surface_model_mpi %s %s %s %s" % (prefix, CS_PATH, gridout, coords, models, output_file)
	print(cmd)
	exitcode = os.system(cmd)
	if (exitcode!=0):
		print("Error generating surface model, aborting.")
		sys.exit(2)
	return output_file

def run_det_smooth_pts(surf_file, coords, nx, ny, smooth_dist):
	'''Usage: %s  <surf file> <model coords file> <nx> <ny> <smoothing dist> <output file>'''
	output_file = "smoothing_pts_mpi.txt"
	cmd = "/bin/date; %s -N 1 -n 1 %s/UCVM/smoothing/determine_smoothing_points.py %s %s %d %d %s %s" % (MPI_CMD, CS_PATH, surf_file, coords, nx, ny, smooth_dist, output_file)
	print(cmd)
	exitcode = os.system(cmd)
	if (exitcode!=0):
		print("Error generating list of smoothing points, aborting.")
		sys.exit(3)
	return output_file

def run_smooth(input_mesh, smooth_pts, nx, ny, nz, smooth_dist, output_mesh):
	'''Usage: %s <AWP-format mesh in> <list of smoothing pts> <RWG nx> <RWG ny> <nz> <smoothing range in pts> <surface_cvm_depth> <velocity mesh out>'''
	if MPI_CMD=="aprun":
                #8 PPN
                np = int(os.environ["PBS_NUM_NODES"])*8
                #Can't run more parallelism than nz
                np = min(nz, np)
                prefix = "%s -n %d -N 8 -S 4" % (MPI_CMD, np)
	elif MPI_CMD=="jsrun":
		np = int(os.environ["LSB_DJOB_NUMPROC"])-1
		nnodes = np/42
		#Can't run more parallelism than nz: divide by 2 because using -a 2
		num_resource_sets = min(nz/2, np)
		#We leave out -r because then different numbers of ranks can be assigned to each processor
		prefix = "%s -a 2 -c 1 -n %d" % (MPI_CMD, num_resource_sets)
	elif MPI_CMD=='srun':
		num_nodes = int(os.environ["SLURM_JOB_NUM_NODES"])
		proc_per_node = int(os.environ["SLURM_CPUS_ON_NODE"])
		#Can't run more parallelism than nz
		np = min(nz, num_nodes * proc_per_node)
		prefix = "%s -N %d -n %d -c 1" % (MPI_CMD, num_nodes, np)
	cmd = "/bin/date; %s %s/UCVM/smoothing/smooth_mpi %s %s %d %d %d %s %s" % (prefix, CS_PATH, input_mesh, smooth_pts, nx, ny, nz, smooth_dist, output_mesh)
	print(cmd)
	exitcode = os.system(cmd)
	if (exitcode!=0):
		print("Error performing smoothing, aborting.")
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
	print("Files gridout, model_coords, input velocity mesh, and output velocity mesh must be specified, aborting.")
	sys.exit(1)

if (smoothing_dist==None):
	print("Smoothing distance in mesh points must be provided, aborting.")
	sys.exit(1)

if (models==None):
	print("Comma-separated list of velocity models must be provided, aborting.")
	sys.exit(1)

#Determine nx, ny, nz from gridout file
with open(gridout, "r") as fp_in:
	data = fp_in.readlines()
	fp_in.close()
	nx = int((data[1].split("="))[1])
	ny = int((data[1+nx+2].split("="))[1])
	nz = int((data[1+nx+2+ny+2].split("="))[1])

#Set LD_LIBRARY_PATH to pick up UCVM libraries
CS_PATH = config.getProperty("CS_PATH")
UCVM_HOME = "%s/UCVM/ucvm_22.7.0_withSFCVM" % CS_PATH
os.environ["LD_LIBRARY_PATH"] = "%s/lib/euclid3/lib:%s/lib/proj/lib:%s/model/cvmsi/lib:%s/model/sfcvm/lib:%s/model/cvms5/lib:%s/model/cca/lib:%s/model/cencal/lib:%s" % (UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME, os.environ["LD_LIBRARY_PATH"])
print(os.environ["LD_LIBRARY_PATH"])

start = default_timer()
surf_model_file = run_det_surf_model(gridout, coords, models)
det_surf_time = default_timer()
smooth_pts_file = run_det_smooth_pts(surf_model_file, coords, nx, ny, smoothing_dist)
det_smooth_pts_time = default_timer()
run_smooth(mesh_in, smooth_pts_file, nx, ny, nz, smoothing_dist, mesh_out)
run_smooth_time = default_timer()

print("Det surf model took %f secs." % (det_surf_time-start))
print("Det smooth pts took %f secs." % (det_smooth_pts_time - det_surf_time))
print("Run smooth took %f secs." % (run_smooth_time - det_smooth_pts_time))
