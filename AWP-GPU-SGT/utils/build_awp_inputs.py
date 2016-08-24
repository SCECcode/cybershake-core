#!/usr/bin/env python

'''This code constructs the 4 inputs needed for running an AWP-ODC-SGT simulation:

1) IN3D file
2) Source description with the source location added
3) Cordfile in AWP format
4) Velocity mesh in AWP format
'''
import sys
import os
import errno

LFS_PATH = "/sw/user/scripts/lfs"

print "Adding config to path."
sys.stdout.flush()
#Add back 3 levels
full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config
import optparse

from build_IN3D import build_IN3D
from build_src import build_src
from build_media import build_media
from build_cordfile import build_cordfile

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as ex:
		if ex.errno==errno.EEXIST and os.path.isdir(path):
			pass
		else: raise


parser = optparse.OptionParser()
parser.add_option("--site", dest="site", action="store", help="Site name")
parser.add_option("--gridout", dest="gridout", action="store", help="Path to gridout input file")
parser.add_option("--fdloc", dest="fdloc", action="store", help="Path to fdloc input file")
parser.add_option("--cordfile", dest="cordfile", action="store", help="Path to cordfile input file")
parser.add_option("--velocity-prefix", dest="vel_prefix", action="store", help="RWG velocity prefix.  If omitted, will not reformat velocity file, just symlink.")
parser.add_option("--frequency", dest="frequency", type=float, action="store", default=0.5, help="Frequency of SGT run, 0.5 Hz by default.")
parser.add_option("--px", dest="px", type=int, action="store", help="Number of processors in X-direction.")
parser.add_option("--py", dest="py", type=int, action="store", help="Number of processors in Y-direction.")
parser.add_option("--pz", dest="pz", type=int, action="store", help="Number of processors in Z-direction.")
parser.add_option("--source-frequency", type=float, dest="source_freq", action="store", help="Low-pass filter frequency to use on the source, default is same frequency as the frequency of the run.")
parser.add_option("--spacing", type=float, dest="spacing", action="store", default=None, help="Override default spacing, derived from frequency.")

(option, args) = parser.parse_args()

site = option.site
gridout = option.gridout
fdloc = option.fdloc
cordfile = option.cordfile

if site==None or gridout==None or fdloc==None or cordfile==None:
	print "site, gridout, fdloc, and cordfile must be specified."
	parser.print_help()
	sys.exit(1)

procs = [option.px, option.py, option.pz]

if (procs[0]==None or procs[1]==None or procs[2]==None):
	print "px, py, pz must be specified."
	parser.print_help()
	sys.exit(1)

rwg_vel_prefix = option.vel_prefix
frequency = option.frequency

source_frequency = frequency
if option.source_freq!=None:
	source_frequency = option.source_freq

awp_comps = ['x', 'y']

for c in awp_comps:
	mkdir_p("comp_%s/input" % c)
	mkdir_p("comp_%s/output_ckp" % c)
	mkdir_p("comp_%s/output_sfc" % c)
	mkdir_p("comp_%s/output_vlm" % c)
	mkdir_p("comp_%s/output_sgt" % c)
	#Set striping for output directory
	exitcode = os.system("%s setstripe -c 160 -s 5m comp_%s/output_sgt" % (LFS_PATH, c))
	if exitcode!=0:
		print "Error striping with command %s setstripe -c 160 -s 5m comp_%s/output_sgt, exiting." % (LFS_PATH, c)
		sys.exit(exitcode)

	print "Building IN3D file for comp %s." % c
	sys.stdout.flush()
	rc = build_IN3D(site, gridout, c, frequency, procs, spacing=option.spacing)
	if not rc==0:
		print "Error in build_IN3D, aborting."
		sys.exit(2)
	print "Building source file."
	sys.stdout.flush()
	rc = build_src(site, fdloc, c, frequency, filter=source_frequency, spacing=option.spacing)
	if not rc==0:
	        print "Error in build_src, aborting."
	        sys.exit(3)
awp_cordfile = "awp.%s.cordfile" % site
print "Building cordfile."
sys.stdout.flush()
#Determine max depth index from gridout file
with open(gridout, "r") as fp_in:
	lines = fp_in.readlines()
	max_depth_index = int(lines[-1].split()[0].strip()) + 1
	fp_in.close()

rc = build_cordfile(site, cordfile, awp_cordfile, max_depth_index)
if not rc==0:
        print "Error in build_cordfile, aborting."
        sys.exit(2)
for c in awp_comps:
	if os.path.lexists("comp_%s/input/%s" % (c, awp_cordfile)):
		os.remove("comp_%s/input/%s" % (c, awp_cordfile))
	os.symlink("../../%s" % awp_cordfile, "comp_%s/input/%s" % (c, awp_cordfile))

awp_media = "awp.%s.media" % (site)
if rwg_vel_prefix is not None:
	print "Building media file."
	sys.stdout.flush()
	rc = build_media(site, gridout, rwg_vel_prefix, awp_media)
	if not rc==0:
	        print "Error in build_media, aborting."
	        sys.exit(2)
else:
	if not os.path.exists(awp_media):
		print "Error, since expected velocity file %s does not exist.  Aborting." % awp_media
		sys.exit(3)
	print "No velocity prefix specified, skipping velocity file reformat."

for c in awp_comps:
	if os.path.lexists("comp_%s/input/%s" % (c, awp_media)):
                os.remove("comp_%s/input/%s" % (c, awp_media))
	os.symlink("../../%s" % awp_media, "comp_%s/input/%s" % (c, awp_media))

