#!/usr/bin/env python

'''This code takes an AWP-ODC-SGT file, reformats it in preparation for post-processing, and constructs the header file needed for running post-processing.'''

import sys
import os

#Add back 3 levels
full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

CS_PATH = config.getProperty("CS_PATH")

def reformat(input_filename, timesteps, output_filename, comp):
	command = "%s/SgtHead/bin/reformat_awp %s %d %s" % (CS_PATH, input_filename, timesteps, output_filename)
	if comp=="z":
		#Add z flag to apply additional factor of two
		command = "%s -z" % command
	print command
	exitcode = os.system(command)
	return exitcode

def write_head(modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, max_f, media, header_name):
	command = "%s/SgtHead/bin/write_head %s %s %s %s %f %d %f %d %s %s %f %s %s" % (CS_PATH, modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, max_f, media, header_name)
	if not comp=="z":
		#Add c flag to use corrected mu
		command = "%s -c " % command
	print command
        exitcode = os.system(command)
        return exitcode


if len(sys.argv)<11:
	print "Usage: %s <site> <AWP SGT> <reformatted SGT filename> <modelbox file> <rwg cordfile> <fdloc file> <gridout file> <IN3D file> <AWP media file> <component> <run_id>" % sys.argv[0]
	sys.exit(0)

SPACING = 0.2
MOMENT = "1.0e20"
MAX_FREQ = 0.5

site = sys.argv[1]
awp_sgt_filename = sys.argv[2]
awp_reformat_sgt_filename = sys.argv[3]
modelbox = sys.argv[4]
cordfile = sys.argv[5]
fdloc = sys.argv[6]
gridout = sys.argv[7]
IN3D = sys.argv[8]
media = sys.argv[9]
comp = sys.argv[10]
run_id = sys.argv[11]

#determine number of output timesteps
fp_in = open(IN3D, "r")
data = fp_in.readlines()
fp_in.close()

params = dict()

for line in data:
	try:
		[value, key] = line.split()
		print key
		params[key] = value
	except ValueError:
		continue

total_ts = int(float(params["TMAX"])/float(params["DT"])+0.5)
decimation = int(params["NTISKP_SGT"])

#rc = reformat(awp_sgt_filename, total_ts/decimation, awp_reformat_sgt_filename, comp)
#if not rc==0:
#	print "Error in reformatting."
#	sys.exit((rc >> 8) & 0xFF)
	
rwg_comp = "z"
if comp=="x":
	rwg_comp = "y"
elif comp=="y":
	rwg_comp = "x"

header_name = "%s_f%s_%s.sgthead" % (site, rwg_comp, run_id)

rc = write_head(modelbox, cordfile, fdloc, gridout, SPACING, total_ts, float(params["DT"]), decimation, comp, MOMENT, MAX_FREQ, media, header_name)
if not rc==0:
        print "Error in header creation."
        sys.exit((rc >> 8) & 0xFF)

sys.exit(0)
