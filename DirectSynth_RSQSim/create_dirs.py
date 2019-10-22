#!/usr/bin/env python2

import sys
import os

'''
This job takes a file with entries like:
src_id/rup_id/rupture_LFN
and creates the directories and symlinks the LFNs.
LFNs look like e<erf_id>_r<rv_id>_<filename on disk>
e51_rv8_1503_18_event200112.srf
'''

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))
sys.path.append(path_add)

import config

if len(sys.argv)<3:
	print "Usage: %s <directory list file> <ERF ID>" % sys.argv[0]
	sys.exit(1)

dirlist = sys.argv[1]
erf_id = int(sys.argv[2])
RUPTURE_ROOT = config.getProperty("RUPTURE_ROOT")
rupture_base_dir = "%s/Ruptures_erf%d" % (RUPTURE_ROOT, erf_id)
with open(dirlist, "r") as fp_in:
	data = fp_in.readlines()
	for line in data:
		pieces = line.split("/")
		src_dir = pieces[0]
		lfn = pieces[1].strip()
		if not os.path.exists(src_dir):
			os.mkdir("%s" % (src_dir), 0755)
		#Strip the erf and rvid from the lfn and put back together
		lfn_pieces = lfn.split("_")
		rup_dir = lfn_pieces[3]
		file_basename = "_".join(lfn_pieces[2:])
		#Figure out where on disk the symlink should come from
		symlink_src = "%s/%s/%s/%s" % (rupture_base_dir, src_dir, rup_dir, file_basename)
		symlink_dest = "%s/%s" % (src_dir, lfn)
		if not os.path.exists(symlink_dest):
			os.symlink(symlink_src, symlink_dest)
	fp_in.close()

