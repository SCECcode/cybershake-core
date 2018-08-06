#!/usr/bin/env python2

'''This script collects Vs30 from UCVM, Vs30 from the velocity mesh, Z1.0, and Z2.5, and populates the database with them.'''

import sys
import os
import MySQLdb
import argparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

parser = argparse.ArgumentParser(description = "'This script collects Vs30 from UCVM, Vs30 from the velocity mesh, Z1.0, and Z2.5, and populates the database with the data.")
parser.add_argument("-lat", "--latitude", help="Latitude of site", type=float)
parser.add_argument("-lon", "--longitude", help="Longitude of site", type=float)
parser.add_argument("-m", "--models", help="Comma-separated string of UCVM velocity models",)
parser.add_argument("-vmesh", "--velocity_mesh", help="Path to velocity mesh file")
parser.add_argument("-fd", "--fdloc", help="Path to fdloc file with site coordinates")
parser.add_argument("-nx", help="Number of grid points in the X direction", type=int)
parser.add_argument("-ny", help="Number of grid points in the Y direction", type=int)
parser.add_argument("-nz", help="Number of grid points in the Z direction", type=int)
parser.add_argument("-s", "--server", help="Database server to use")
parser.add_argument("-r", "--run_id", help="Run ID")

args = parser.parse_args()
if args.latitude==None or args.longitude==None or args.models==None or args.velocity_mesh==None or args.nx==None or args.ny==None or args.nz==None or args.fdloc==None or args.run_id==None:
	print "Missing required arguments."
	parser.print_help()
	sys.exit(1)

server = "moment.usc.edu"
if args.server is not None:
	server = args.server

#Call get_model_info_for_db to get model Vs30, Z1.0, Z2.5
model_output_file = "ucvm_data.txt"
cmd = "%s/UCVM/bin/get_model_info_for_db %f %f %s %s" % (config.getProperty("CS_PATH")), args.latitude, args.longitude, args.models, model_output_file)
os.system(cmd)

#Run value_from_mesh
with open(args.fdloc, "r") as fp_in:
	[x_index, y_index] = fp_in.readline().split()
	fp_in.close()
z_index = 0
mesh_output_file = "mesh_data.txt"
cmd = "%s/UCVM/utils/value_from_mesh %s %d %d %d %d %d %d > %s" % (args.velocity_mesh, args.nx, args.ny, args.nz, x_index, y_index, z_index, mesh_output_file)
os.system(cmd)

with open(model_output_file, "r") as fp_in:
	model_vs30 = float(fp_in.readline())
	model_z10 = float(fp_in.readline())
	model_z25 = float(fp_in.readline())
	fp_in.close()

with open(mesh_output_file, "r") as fp_in:
	[mesh_vp, mesh_vs, mesh_rho] = [float(i) for i in fp_in.readline().split()]
	fp_in.close()

conn = MySQLdb.connect(host=server, db="CyberShake", user="cybershk", passwd='***REMOVED***')
cur = conn.cursor()
update = "update CyberShake_Runs set Model_Vs30=%f, Mesh_Vs30=%f, Z1_0=%f, Z2_5=%f where Run_ID=%d" % (model_vs30, mesh_vs30, model_z10, model_z25, run_id)
cur.query(update)

