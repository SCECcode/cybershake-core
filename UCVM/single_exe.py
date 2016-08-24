#!/usr/bin/env python

import os
import sys
import optparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

parser = optparse.OptionParser()
parser.add_option("--site", dest="site", action="store", help="Site name")
parser.add_option("--gridout", dest="gridout", action="store", help="Path to gridout (output)")
parser.add_option("--coordfile", dest="coordsfile", action="store", help="Path to coordfile (output)")
parser.add_option("--models", dest="models", action="store", help="Comma-separated string on velocity models to use.")
parser.add_option("--format", dest="format", action="store", help="Specify awp or rwg format for output.")
parser.add_option("--frequency", dest="frequency", type="float", action="store", help="Frequency")
parser.add_option("--spacing", dest="spacing", type="float", action="store", help="Override default spacing with this value")
parser.add_option("--min_vs", dest="min_vs", type="float", action="store", help="Override minimum Vs value.  Minimum Vp and minimum density will be 3.4 times this value.")

(option, args) = parser.parse_args()

site = option.site
gridout = option.gridout
modelcords = option.coordsfile
models = option.models
format = option.format
if site==None or gridout==None or modelcords==None or models==None or format==None:
	print "All of site, gridout, modelcords, models, and format must be specified."
	parser.print_help()
	sys.exit(-1)

if option.frequency is not None:
	frequency = option.frequency
else:
	frequency = 0.5

if option.spacing is not None:
	zstep = option.spacing
else:
	zstep = 100.0/frequency

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

#set up for striping if awp
if format=="awp":
	os.system("/usr/bin/lfs setstripe -c 100 -s 5m awp.%s.media" % site)

if option.min_vs is not None:
        command = '%s/single_exe.csh %s %s %d %d %d %s %s %s %s %s %s %s %f' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id, format, option.min_vs)
else:
	command = '%s/single_exe.csh %s %s %d %d %d %s %s %s %s %s %s %s' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id, format)
print command
exitcode = os.system(command)
sys.exit((exitcode >> 8) & 0xFF)
