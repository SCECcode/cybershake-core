#!/usr/bin/env python3

import sys
import os

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))
sys.path.append(path_add)

import config

MPI_CMD = config.getProperty("MPI_CMD")

if len(sys.argv)<2:
	print("Usage: %s <file>")
	sys.exit(1)

filename = sys.argv[1]

if MPI_CMD=="aprun":
	cmd = "aprun -n 1 md5sum %s > %s.md5" % (filename, filename)
elif MPI_CMD=="jsrun":
	cmd = "jsrun -n 1 /usr/bin/md5sum %s > %s.md5" % (filename, filename)
else:
	print("Don't recognize mpi cmd %s, aborting." % MPI_CMD)
	sys.exit(1)

os.system(cmd)

