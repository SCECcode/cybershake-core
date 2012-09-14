#!/usr/bin/env python

import os
import sys

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

if len(sys.argv)<5:
	print 'Usage: sim_sgt.py <site> <gridout> <modelbox> <fdloc> <component>'
	sys.exit(0)

site = sys.argv[1]
gridout = sys.argv[2]
modelbox = sys.argv[3]
fdloc = sys.argv[4]
cordfile = sys.argv[5]
component = sys.argv[6]

#get grid steps
input = open(gridout)
gridoutContents = input.readlines()
input.close()
ns = []
ns.append(int((gridoutContents[1].split("="))[1]))
ns.append(int((gridoutContents[1+ns[0]+2].split("="))[1]))
ns.append(int((gridoutContents[1+ns[0]+2+ns[1]+2].split("="))[1]))

#get mlat, mlon, mrot
input = open(modelbox)
modelboxContents = input.readlines()
input.close()
modelTokens = modelboxContents[4].split()
mlon = (float)(modelTokens[1])
mlat = (float)(modelTokens[3])
mrot = (float)(modelTokens[5])

#get fdloc
input = open(fdloc)
fdLine = input.readline()
tokens = fdLine.split()
xsrc = (int)(tokens[0])
ysrc = (int)(tokens[1])
input.close()

#get config
cs_path = config.getProperty('CS_PATH')
scratch_path = config.getProperty('SCRATCH_PATH')
tmp_path = config.getProperty('TMP_PATH')
mpi_cmd = config.getProperty('MPI_CMD')
job_id = config.getJobID()

command = '%s/sim_sgt.csh %s %d %d %d %f %f %f %d %d %s %s %s %s %s %s' % (sys.path[0], site, ns[0], ns[1], ns[2], mlat, mlon, mrot, xsrc, ysrc, component, cs_path, scratch_path, tmp_path, mpi_cmd, job_id)
print command
exitcode = os.system(command)
sys.exit((exitcode >> 8) & 0xFF)

