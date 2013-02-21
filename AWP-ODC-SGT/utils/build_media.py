#!/usr/bin/env python

import sys
import os
import config

def build_media(site, gridout, rwg_vel_prefix, media_out):
	mpi_cmd = config.getProperty("MPI_CMD")

	num_nodes = int(os.getenv("PBS_NUM_NODES"))
	num_ppn = int(os.getenv("PBS_NUM_PPN"))
	np = num_nodes*num_ppn
	
	cs_path = config.getProperty("CS_PATH")
	
	#determine NX, NY, NZ
	fp_in = open(gridout, "r")
	data = fp_in.readlines()
	fp_in.close()
	nx = int((data[1].split("="))[1])
	ny = int((data[1+nx+2].split("="))[1])
	nz = int((data[1+nx+2+ny+2].split("="))[1])
	
	cmd = "%s -n %d %s/UCVM/bin/reformat_velocity_mpi %d %d %d %s %s" % (mpi_cmd, np, cs_path, nx, ny, nz, rwg_vel_prefix, media_out)
	
	print cmd
	exitcode = os.system(cmd)
	return ((exitcode >> 8) & 0xFF)

