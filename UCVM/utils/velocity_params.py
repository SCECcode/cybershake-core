#!/usr/bin/env python

'''This script collects Vs30 from UCVM, Vs30 from the velocity mesh, Z1.0, and Z2.5, and populates the database with them.'''

import sys
import os
import pymysql
import argparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

parser = argparse.ArgumentParser(description = "'This script collects Vs30 from UCVM, Vs30 from the velocity mesh, Z1.0, and Z2.5, and populates the database with the data.")
parser.add_argument("-lat", "--latitude", help="Latitude of site", type=float)
parser.add_argument("-lon", "--longitude", help="Longitude of site", type=float)
parser.add_argument("-m", "--models", help="Comma-separated string of UCVM velocity models",)
parser.add_argument("-vmesh", "--velocity_mesh", help="Path to velocity mesh file")
parser.add_argument("-fd", "--fdloc", help="Path to fdloc file with site coordinates")
parser.add_argument("-go", "--gridout", help="Path to gridout file")
parser.add_argument("-s", "--server", help="Database server to use")
parser.add_argument("-r", "--run_id", help="Run ID")

args = parser.parse_args()
if args.latitude==None or args.longitude==None or args.models==None or args.velocity_mesh==None or args.gridout==None or args.fdloc==None or args.run_id==None:
	print("Missing required arguments.")
	parser.print_help()
	sys.exit(1)

server = "moment.usc.edu"
if args.server is not None:
	server = args.server

#Call get_model_info_for_db to get model Vs30, Z1.0, Z2.5
model_output_file = "ucvm_data.txt"
UCVM_HOME = "/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0"
ld_lib_path="%s/lib/euclid3/lib:%s/lib/proj-4/lib:%s/model/cvms426/lib:%s/model/cencal/lib:%s/model/cvms5/lib:%s/model/cca/lib" % (UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME, UCVM_HOME)
proj_lib_path="/gpfs/alpine/proj-shared/geo112/CyberShake/software/UCVM/ucvm-18.5.0/lib/proj-4/share/proj"
cmd = "export PROJ_LIB=%s; export LD_LIBRARY_PATH=%s:LD_LIBRARY_PATH; %s/UCVM/bin/get_model_info_for_db %f %f %s %s" % (proj_lib_path, ld_lib_path, config.getProperty("CS_PATH"), args.latitude, args.longitude, args.models, model_output_file)
print(cmd)
os.system(cmd)

#Determine mesh size
with open(args.gridout, "r") as fp_in:
	data = fp_in.readlines()
	nx = int(data[1].split("=")[1].strip())
	ny = int(data[1+nx+2].split("=")[1].strip())
	nz = int(data[1+nx+2+ny+2].split("=")[1].strip())
	print(nx, ny, nz)
	fp_in.close()

#Run value_from_mesh
with open(args.fdloc, "r") as fp_in:
	[x_index, y_index] = [int(x) for x in fp_in.readline().split()]
	fp_in.close()
z_index = 0
mesh_output_file = "mesh_data.txt"
cmd = "%s/UCVM/utils/value_from_mesh %s %d %d %d %d %d %d > %s" % (config.getProperty("CS_PATH"), args.velocity_mesh, nx, ny, nz, x_index, y_index, z_index, mesh_output_file)
print(cmd)
os.system(cmd)

with open(model_output_file, "r") as fp_in:
	model_vs30 = float(fp_in.readline())
	model_z10 = float(fp_in.readline())
	model_z25 = float(fp_in.readline())
	fp_in.close()

with open(mesh_output_file, "r") as fp_in:
	[mesh_vp, mesh_vs, mesh_rho] = [float(i) for i in fp_in.readline().split()]
	fp_in.close()

<<<<<<< HEAD
conn = pymysql.connect(host=server, db="CyberShake", user="cybershk", passwd='re@lStil1')
=======
with open(args.password_file, "r") as fp_in:
    pieces = fp_in.readline().split(":")
    username = pieces[0]
    passwd = pieces[1]
    fp_in.close()

conn = pymysql.connect(host=server, db="CyberShake", user=username, passwd=passwd)
>>>>>>> c27aacadfa2ae8283058f6863a97383aaeefc65f
cur = conn.cursor()
update = "update CyberShake_Runs set Model_Vs30=%f, Mesh_Vsitop=%f, Z1_0=%f, Z2_5=%f where Run_ID=%d" % (model_vs30, mesh_vs, model_z10, model_z25, int(args.run_id))
print(update)
cur.execute(update)
conn.commit()
cur.close()
