#!/usr/bin/env python

import sys;
import os;

if len(sys.argv) < 9:
    print "Usage: pre_cvm.py <site> <erf_id> <modelboxFile> <gridfile> <gridout> <coordfile> <paramsfile> <boundsfile> [-frequency <freq] [-gpu]"
    print "Example: pre_cvm.py USC 34 USC.modelbox gridfile_USC gridout_USC model_coords_GC_USC model_params_GC_USC model_bounds_GC_USC"
    sys.exit(-1)

site = sys.argv[1]
erfID = sys.argv[2]
modelbox = os.path.abspath(sys.argv[3])
gridfile = os.path.abspath(sys.argv[4])
gridout = os.path.abspath(sys.argv[5])
coordsfile = os.path.abspath(sys.argv[6])
paramsfile = os.path.abspath(sys.argv[7])
boundsfile = os.path.abspath(sys.argv[8])
frequency = 0.5
gpu_arg = ""
if len(sys.argv) >= 10:
	i = 9
	while i<len(sys.argv):
		if sys.argv[i]=="-gpu":
			#Using GPU rules; adjust volume size
			print "Using GPU volume rules."
			gpu_arg = "gpu"
		elif sys.argv[i]=="-frequency":
			frequency = float(sys.argv[i+1])
			i += 1
		i += 1

os.chdir(os.path.join(sys.path[0], "Modelbox"))
exitcode = os.system("./get_modelbox.py %s %s %s %s" % (site, erfID, modelbox, gpu_arg))
if exitcode!=0:
	sys.exit((exitcode >> 8) & 0xFF)
os.chdir("../GenGrid_py")
exitcode = os.system("./gen_grid.py %s %s %s %s %s %s %f %s" % (modelbox, gridfile, gridout, coordsfile, paramsfile, boundsfile, frequency, gpu_arg))
sys.exit((exitcode >> 8) & 0xFF)
