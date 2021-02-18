#!/usr/bin/env python

import sys
import os

def build_src(site, fdloc, awp_comp, frequency, nt, filter=None, spacing=None):
	if awp_comp=='x':
		comp = "y"
	elif awp_comp=='y':
		comp = "x"
	elif awp_comp=='z':
		comp = "z"
	else:
		print "Error:  component %s not recognized, aborting." % comp
		sys.exit(1)
	
	fp_in = open(fdloc, "r")
	data = fp_in.readline()
	fp_in.close()
	
	[src_x, src_y] = data.split()
	
	#Add to support a filter frequency different from the simulation frequency
	if filter==None:
		filter = frequency

	source_name = "%s/data/f%s_src_%d_%.1fhzFilter" % (sys.path[0], awp_comp, nt, filter)
	#print("Using source %s" % source_name)
	if not os.path.exists(source_name):
		print "Error: could not find source file %s with nt=%d and filter frequency = %.1f, aborting." % (source_name, nt, filter)
		return 1

	fp_in = open("%s" % (source_name), "r")

	#fp_in = open("%s/data/f%s_src_%d_2hzFilter" % (sys.path[0], awp_comp, nt))
	#fp_in = open("%s/data/f%s_src_%d" % (sys.path[0], awp_comp, nt))
	src_data = fp_in.readlines()
	fp_in.close()
	
	fp_out = open("comp_%s/input/%s_f%s_src" % (awp_comp, site, awp_comp), "w")
	#swap X and Y, add 1 to each since AWP uses 1-indexing
	fp_out.write("%d %d 1\n" % (int(src_y)+1, int(src_x)+1))
	for line in src_data:
		fp_out.write(line)
	fp_out.flush()
	fp_out.close()
	return 0

if __name__=="__main__":
	#build_src(site, fdloc, awp_comp, frequency, nt, filter=None, spacing=None)
	if len(sys.argv)<6:
		print("Usage: %s <site> <fdloc file> <awp component> <frequency for AWP> <nt for AWP> [filter frequency]" % sys.argv[0])
		sys.exit(1)

	site=sys.argv[1]
	fdloc_filename=sys.argv[2]
	awp_component=sys.argv[3]
	freq=float(sys.argv[4])
	nt=int(sys.argv[5])
	filter_freq = freq
	if len(sys.argv)>6:
		filter_freq=float(sys.argv[6])
	sys.exit(build_src(site, fdloc_filename, awp_component, freq, nt, filter=filter_freq, spacing=None))


