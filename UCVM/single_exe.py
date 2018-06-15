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
parser.add_option("--spacing", dest="spacing", type="float", action="store", help="Override default spacing with this value (km)")
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
	zstep = option.spacing*1000.0
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

#if cca is one of the models, check to see if GTL is on
with open('%s/UCVM/ucvm-18.5.0/model/cca/data/config' % cs_path, 'r') as fp_in:
	data = fp_in.readlines()
	for line in data:
		if line.find('gtl')>-1:
			pieces = line.split('=')
			if pieces[1].strip()=='off':
				print "CCA: GTL is off."
			elif pieces[1].strip()=='on':
				print "CCA: GTL is ON."
			else:
				print "CCA: GTL status unknown."
			break
	fp_in.close()

#set up for striping if awp
LFS_CMD = "/opt/cray/lustre-cray_gem_s/2.5_3.0.101_0.31.1_1.0502.8394.10.1-1.0502.17198.8.50/bin/lfs"
if format=="awp":
	os.system("%s setstripe -c 100 -s 5m awp.%s.media" % (LFS_CMD, site))

if option.min_vs is not None:
        command = '%s/single_exe.csh %s %s %d %d %d %s %s %s %s %s %s %s %.01f' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id, format, option.min_vs)
else:
	command = '%s/single_exe.csh %s %s %d %d %d %s %s %s %s %s %s %s' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id, format)
print command
exitcode = os.system(command)
sys.exit((exitcode >> 8) & 0xFF)
