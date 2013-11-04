#!/usr/bin/env python

import os
import sys

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

if len(sys.argv)<5:
	print 'Usage: genmod.py <site> <gridout> <model_cords file> <models> [frequency]'
	sys.exit(-1)

site = sys.argv[1]
gridout = sys.argv[2]
modelcords = sys.argv[3]
models = sys.argv[4]
if len(sys.argv)==6:
	frequency = float(sys.argv[5])
	zstep = 100.0/frequency #meters

#get grid steps
input = open(gridout)
gridoutContents = input.readlines()
input.close()
ns = []
ns.append(int((gridoutContents[1].split("="))[1]))
ns.append(int((gridoutContents[1+ns[0]+2].split("="))[1]))
ns.append(int((gridoutContents[1+ns[0]+2+ns[1]+2].split("="))[1]))

#generate zdepths file
zdepths = open('zdepths.meter', 'w')
for i in range(ns[2]):
	zdepths.write('%.2f\n' % (i * zstep))

zdepths.flush()
zdepths.close()

cs_path = config.getProperty("CS_PATH")
scratch_path = config.getProperty("SCRATCH_PATH")
log_root = config.getProperty("LOG_PATH")
mpi_cmd = config.getProperty("MPI_CMD")
job_id = config.getJobID()

command = '%s/genmod.csh %s %s %d %d %d %s %s %s %s %s %s' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id)
print command
exitcode = os.system(command)
sys.exit((exitcode >> 8) & 0xFF)
