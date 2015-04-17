#!/usr/bin/env python

import sys
import os

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))
sys.path.append(path_add)

import config

#Start memcached on each node
memcached_path = "/projects/sciteam/jmz/CyberShake/utils/pegasus_wrappers/invoke_memcached.sh"
cmd = "aprun -N 1 %s" % memcached_path
rc = os.system(cmd)
if rc!=0:
	print "Error launching memcached.  Will continue with DirectSynth anyway."

#Launch direct_synth
ds_path = "%s/bin/direct_synth %s" % (sys.path[0], " ".join(sys.argv[1:]))
num_nodes = int(os.environ['PBS_NUM_NODES'])
ppn = int(os.environ['PBS_NUM_PPN'])
cmd = "aprun -n %d -N %d %s" % (num_nodes*ppn, ppn, ds_path)
print cmd
rc = os.system(cmd)
sys.exit((rc >> 8) & 0xFF)
