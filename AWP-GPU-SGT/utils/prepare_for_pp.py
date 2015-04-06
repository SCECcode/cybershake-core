#!/usr/bin/env python

'''This code takes an AWP-ODC-SGT file, reformats it in preparation for post-processing, and constructs the header file needed for running post-processing.'''

import sys
import os
import optparse

#Add back 3 levels
full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

CS_PATH = config.getProperty("CS_PATH")

def reformat(input_filename, timesteps, num_pts, output_filename, comp):
	#On Titan, need aprun to execute this on a compute node, not the aprun node
	command = "aprun -n 4 -N 2 -S 1 %s/SgtHead/bin/reformat_awp_mpi %s %d %d %s" % (CS_PATH, input_filename, timesteps, num_pts, output_filename)
	if comp=="z":
		#Add z flag to apply additional factor of two
		command = "%s -z" % command
	print command
	sys.stdout.flush()
	exitcode = os.system(command)
	return exitcode

def write_head(modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, source_freq, media, header_name):
	#On Titan, need aprun to execute this on a compute node, not the aprun node
	command = "aprun -n 1 %s/SgtHead/bin/write_head %s %s %s %s %f %d %f %d %s %s %f %s %s" % (CS_PATH, modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, source_freq, media, header_name)
	if not comp=="z":
		#Add c flag to use corrected mu
		command = "%s -c " % command
	print command
	sys.stdout.flush()
        exitcode = os.system(command)
        return exitcode

usage = "Usage: %s <site> <AWP SGT> <reformatted SGT filename> <modelbox file> <rwg cordfile> <fdloc file> <gridout file> <IN3D file> <AWP media file> <component> <run_id> <header> [frequency]" % sys.argv[0]

if len(sys.argv)<13:
	print usage
	sys.exit(1)

MOMENT = "1.0e20"

parser = optparse.OptionParser(usage = usage)
parser.add_option("-n", "--no-md5", dest="no_md5", action="store_true", default=False, help="Skip MD5 sum step.")
parser.add_option("-s", "--source-frequency", type="float", dest="src_freq", action="store", help="Frequency of the SGT source (default is same as simulation frequency)")
(options, args) = parser.parse_args()
skip_md5 = options.no_md5
source_freq = options.src_freq

site = args[0]
awp_sgt_filename = args[1]
awp_reformat_sgt_filename = args[2]
modelbox = args[3]
cordfile = args[4]
fdloc = args[5]
gridout = args[6]
IN3D = args[7]
media = args[8]
comp = args[9]
run_id = args[10]
header_out = args[11]
MAX_FREQ = float(args[12])
print "Max frequency: %f" % (MAX_FREQ)

if source_freq==None:
	source_freq = MAX_FREQ

#determine number of output timesteps
fp_in = open(IN3D, "r")
data = fp_in.readlines()
fp_in.close()

#determine number of points
fp_in = open(cordfile, "r")
fp_in.readline()
fp_in.readline()
fp_in.readline()
fp_in.readline()
num_sgt_pts = int(fp_in.readline())
fp_in.close()

params = dict()

for line in data:
	try:
		[value, key] = line.split()
		params[key] = value
	except ValueError:
		continue

total_ts = int(float(params["TMAX"])/float(params["DT"])+0.5)
decimation = int(params["NTISKP_SGT"])
spacing = float(params["DH"])/1000.0

rc = reformat(awp_sgt_filename, total_ts/decimation, num_sgt_pts, awp_reformat_sgt_filename, comp)
rc = 0
if not rc==0:
	print "Error in reformatting."
	sys.exit((rc >> 8) & 0xFF)
	
rwg_comp = "z"
if comp=="x":
	rwg_comp = "y"
elif comp=="y":
	rwg_comp = "x"

header_name = "%s_f%s_%s.sgthead" % (site, rwg_comp, run_id)

rc = write_head(modelbox, cordfile, fdloc, gridout, spacing, total_ts, float(params["DT"]), decimation, rwg_comp, MOMENT, source_freq, media, header_name)
if not rc==0:
        print "Error in header creation."
        sys.exit((rc >> 8) & 0xFF)

#Run md5sum
if not skip_md5:
	os.system("md5sum %s > %s.md5" % (awp_reformat_sgt_filename, awp_reformat_sgt_filename))

sys.exit(0)
