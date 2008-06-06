#!/u/ac/scottcal/python-2.5.1/bin/python

import sys;
import os;

if len(sys.argv) < 8:
    print "Usage: pre_cvm.py <site> <modelboxFile> <gridfile> <gridout> <coordfile> <paramsfile> <boundsfile>"
    print "Example: pre_cvm.py USC USC.modelbox gridfile_USC gridout_USC model_coords_GC_USC model_params_GC_USC model_bounds_GC_USC"
    sys.exit(-1)

site = sys.argv[1]
modelbox = os.path.abspath(sys.argv[2])
gridfile = os.path.abspath(sys.argv[3])
gridout = os.path.abspath(sys.argv[4])
coordsfile = os.path.abspath(sys.argv[5])
paramsfile = os.path.abspath(sys.argv[6])
boundsfile = os.path.abspath(sys.argv[7])

os.chdir(os.path.join(sys.path[0], "Modelbox"))
exitcode = os.system("./get_modelbox.py %s %s" % (site, modelbox))
if exitcode!=0:
	sys.exit((exitcode >> 8) & 0xFF)
os.chdir("../GenGrid_py")
exitcode = os.system("./gen_grid.py %s %s %s %s %s %s" % (modelbox, gridfile, gridout, coordsfile, paramsfile, boundsfile))
sys.exit((exitcode >> 8) & 0xFF)
