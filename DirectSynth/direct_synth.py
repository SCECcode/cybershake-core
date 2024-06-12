#!/usr/bin/env python

import sys
import os

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))
sys.path.append(path_add)

import config

MPI_CMD = config.getProperty("MPI_CMD")

if MPI_CMD=="aprun":
	num_nodes = int(os.environ['PBS_NUM_NODES'])
	ppn = int(os.environ['PBS_NUM_PPN'])
	#Start memcached on each node
	memcached_path = "/projects/sciteam/bahm/CyberShake/utils/pegasus_wrappers/invoke_memcached.sh"
	cmd = "aprun -N 1 -n %d %s" % (num_nodes, memcached_path)
elif MPI_CMD=="jsrun":
	num_res_sets = int(os.environ['LSB_DJOB_NUMPROC'])-1
	num_nodes = num_res_sets/42
	#1 is taken for memcached
	#Start memcached once on each node
	memcached_path = "/gpfs/alpine/proj-shared/geo112/CyberShake/utils/pegasus_wrappers/invoke_memcached.sh"
	cmd = "jsrun -a 1 -c 1 -g 0 -r 1 -n %d %s" % (num_nodes, memcached_path)
else:
	print("Don't know what to do with MPI_CMD=%s, aborting." % MPI_CMD)
	sys.exit(1)
print cmd
rc = os.system(cmd)
if rc!=0:
        print "Error launching memcached.  Will continue with DirectSynth anyway."

#Launch direct_synth
ds_path = "%s/bin/direct_synth %s" % (sys.path[0], " ".join(sys.argv[1:]))
if MPI_CMD=="aprun":
	cmd = "export MPICH_MPIIO_HINTS=*:romio_cb_read=disable; export APRUN_BALANCED_INJECTION=64; ulimit -c unlimited; export ATP_ENABLED=1; aprun -n %d -N %d %s" % (num_nodes*ppn, ppn, ds_path)
elif MPI_CMD=="jsrun":
	#Can't have more than 42 resource sets/core
	#cmd = "jsrun -n %d -a 4 -c 1 -g 0 -r 42 %s" % (num_res_sets, ds_path)
	cmd = "jsrun -n %d -a 1 -c 1 -g 0 -r 42 %s" % (num_res_sets, ds_path)
	cmd = "which jsrun; module list; %s" % (cmd)
elif MPI_CMD=="ibrun":
	#cmd = "export IBRUN_TASKS_PER_NODE=56; ibrun %s" % ds_path
    cmd = "export LD_LIBRARY_PATH=/work2/00349/scottcal/frontera/CyberShake/utils/libmemcached_1.0.18/lib:$LD_LIBRARY_PATH; ibrun %s" % ds_path
print(cmd)
rc = os.system(cmd)
sys.exit((rc >> 8) & 0xFF)

