#!/usr/bin/env python

import os
import sys

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

if len(sys.argv)<5:
	print 'Usage: merge_cvm.py <site> <gridfile> <velocity p-file> <velocity s-file> <velocity d-file>'
	sys.exit(-1)

site = sys.argv[1]
gridfile = sys.argv[2]
p_file = sys.argv[3]
s_file = sys.argv[4]
d_file = sys.argv[5]

cs_path = config.getProperty("CS_PATH")
scratch_path = config.getProperty("SCRATCH_PATH")
log_root = config.getProperty("LOG_PATH")
mpi_cmd = config.getProperty("MPI_CMD")
job_id  = config.getJobID()

command = '%s/merge_cvm.csh %s %s %s %s %s %s %s %s %s %s' % (sys.path[0], site, gridfile, p_file, s_file, d_file, cs_path, scratch_path, log_root, mpi_cmd, job_id)
print command
exitcode = os.system(command)
sys.exit((exitcode >> 8) & 0xFF)
