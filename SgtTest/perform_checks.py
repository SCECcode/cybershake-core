#!/usr/bin/env python3

'''This performs a check of SGT size, then calls the nan check code.'''

import sys
import os
import struct

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

if len(sys.argv)<4:
	print("Usage: %s <SGT file> <cordfile> <IN3D file>" % sys.argv[0])
	sys.exit(1)

sgt_filename = sys.argv[1]
cordfile_name = sys.argv[2]
in3d_filename = sys.argv[3]

#Get actual file size
sgt_filesize = os.path.getsize(sgt_filename)

#Read out header information in header file to check for match

with open(cordfile_name, "r") as fp_in:
	line = fp_in.readline()
	while line[0]=="#":
		line = fp_in.readline()
	np = int(line.strip())
	fp_in.close()

#Determine number of timesteps and timeskip
with open(in3d_filename, "r") as fp_in:
	data = fp_in.readlines()
	for line in data:
		pieces = line.split()
		if len(pieces)>1 and pieces[1].strip()=="NST":
			nt = int(pieces[0].strip())
			continue
		if len(pieces)>1 and pieces[1].strip()=="NTISKP_SGT":
			timeskip = int(pieces[0].strip())
			continue
	fp_in.close()

FLOAT_SIZE = 4
SGT_COMPONENTS = 6

#Divide nt by timeskip, because we decimate in time
expected_size = np*nt/timeskip*SGT_COMPONENTS*FLOAT_SIZE

if expected_size!=sgt_filesize:
	print("Error: cordfile %s leads us to expect %d points and IN3D file %s has nt=%d, for an SGT filesize of %d, but the SGT file actually has size %d.  Aborting." % (cordfile_name, np, in3d_filename, nt, expected_size, sgt_filesize))
	sys.exit(2)

cs_path = config.getProperty("CS_PATH")
mpi_cmd = config.getProperty("MPI_CMD")

#Make sure it gets executed on a compute node
if mpi_cmd=="aprun":
	mpi_cmd = "aprun -n 1"
elif mpi_cmd=="mpirun":
	mpi_cmd = "mpirun -np 1"

cmd = "%s %s/SgtTest/bin/check_for_nans %s" % (mpi_cmd, cs_path, sgt_filename)
exitcode = os.system(cmd)
if exitcode!=0:
	print("Error when checking for NaNs, zeros.")
	sys.exit(3)
sys.exit(exitcode)
