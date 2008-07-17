#!/usr/bin/python

import sys;
import os;

if len(sys.argv) < 9:
    print "Usage: pre_cvm.py <site> <erf_id> <modelboxFile> <gridfile> <gridout> <coordfile> <paramsfile> <boundsfile>"
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

os.chdir(os.path.join(sys.path[0], "Modelbox"))
exitcode = os.system("./get_modelbox.py %s %s %s" % (site, erfID, modelbox))
if exitcode!=0:
	sys.exit((exitcode >> 8) & 0xFF)
os.chdir("../GenGrid_py")
exitcode = os.system("./gen_grid.py %s %s %s %s %s %s" % (modelbox, gridfile, gridout, coordsfile, paramsfile, boundsfile))
sys.exit((exitcode >> 8) & 0xFF)
