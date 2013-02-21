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

#Add back 3 levels
full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

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


if len(sys.argv)<6:
	print "Usage: %s <site> <gridout> <rwg velocity prefix> <fdloc> <rwg cordfile>" % (sys.argv[0])
	sys.exit(1)

site = sys.argv[1]
gridout = sys.argv[2]
rwg_vel_prefix = sys.argv[3]
fdloc = sys.argv[4]
rwg_cordfile = sys.argv[5]

awp_comps = ['x', 'y']

for c in awp_comps:
	mkdir_p("comp_%s/input" % c)
	mkdir_p("comp_%s/output_ckp" % c)
	mkdir_p("comp_%s/output_sfc" % c)
	mkdir_p("comp_%s/output_vlm" % c)
	mkdir_p("comp_%s/output_sgt" % c)

	rc = build_IN3D(site, gridout, c)
	if not rc==0:
		print "Error in build_IN3D, aborting."
		sys.exit(2)
	rc = build_src(site, fdloc, c)
	if not rc==0:
	        print "Error in build_src, aborting."
	        sys.exit(3)
awp_cordfile = "awp.%s.cordfile" % site
rc = build_cordfile(site, rwg_cordfile, awp_cordfile)
if not rc==0:
        print "Error in build_cordfile, aborting."
        sys.exit(2)
for c in awp_comps:
	if os.path.exists("comp_%s/input/%s" % (c, awp_cordfile)):
		os.remove("comp_%s/input/%s" % (c, awp_cordfile))
	os.symlink("../../%s" % awp_cordfile, "comp_%s/input/%s" % (c, awp_cordfile))
awp_media = "awp.%s.media" % (site)
rc = build_media(site, gridout, rwg_vel_prefix, awp_media)
if not rc==0:
        print "Error in build_media, aborting."
        sys.exit(2)
for c in awp_comps:
	if os.path.exists("comp_%s/input/%s" % (c, awp_media)):
                os.remove("comp_%s/input/%s" % (c, awp_media))
	os.symlink("../../%s" % awp_media, "comp_%s/input/%s" % (c, awp_media))

