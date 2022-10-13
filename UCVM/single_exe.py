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
parser.add_option("--h_fraction", dest="h_frac", type="float", action="store", help="Depth, in fractions of a grid point, to query UCVM at when populating the surface points.")
parser.add_option("--ely-taper", dest="ely_taper", type="str", action="store", help="Type of Ely taper to use.  Choices are 'none' (the default), 'all' (apply the taper to all points), or 'ifless' (apply the taper at all points for which the taper has a lower Vs)")
parser.add_option("--taper-depth", dest="taper_depth", type="float", action="store", help="Transition depth in meters for the Ely taper, if it's being applied.  Default is 700m.")

(option, args) = parser.parse_args()

site = option.site
gridout = option.gridout
modelcords = option.coordsfile
models = option.models
format = option.format
if site==None or gridout==None or modelcords==None or models==None or format==None:
	print("All of site, gridout, modelcords, models, and format must be specified.")
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

if option.ely_taper is not None:
    ely_taper = option.ely_taper
else:
    ely_taper = "none"

if option.taper_depth is not None:
    taper_depth = option.taper_depth
else:
    taper_depth = 700.0


#Depth to query UCVM at for surface points, in m
surface_cvm_depth = 0.0
if option.h_frac is not None:
	surface_cvm_depth = float(option.h_frac)*zstep

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
with open('%s/UCVM/ucvm-22.7.0/model/cca/data/config' % cs_path, 'r') as fp_in:
	data = fp_in.readlines()
	for line in data:
		if line.find('gtl')>-1:
			pieces = line.split('=')
			if pieces[1].strip()=='off':
				print("CCA: GTL is off.")
			elif pieces[1].strip()=='on':
				print("CCA: GTL is ON.")
			else:
				print("CCA: GTL status unknown.")
			break
	fp_in.close()

if option.min_vs is not None:
        command = '%s/single_exe.csh %s %s %d %d %d %s %s %s %s %s %s %s %.1f %s %f %.01f' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id, format, surface_cvm_depth, ely_taper, taper_depth, option.min_vs)
else:
	command = '%s/single_exe.csh %s %s %d %d %d %s %s %s %s %s %s %s %.1f %s %f' % (sys.path[0], site, modelcords, ns[0], ns[1], ns[2], cs_path, scratch_path, log_root, models, mpi_cmd, job_id, format, surface_cvm_depth, ely_taper, taper_depth)
print(command)
exitcode = os.system(command)
sys.exit((exitcode >> 8) & 0xFF)
