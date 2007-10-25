#!/usr/bin/python

import sys;
import os;

if len(sys.argv) < 2:
    print "Usage: pre_cvm.py <site>"
    print "Example: pre_cvm.py USC"
    sys.exit()

site = sys.argv[1]
modelboxFile = site + ".modelbox"
os.system("mkdir ../../data/ModelParams/" + site)
os.chdir("Modelbox")
os.system("./get_modelbox.py " + site + " ../../../data/ModelParams/" + site + "/" + modelboxFile)
os.chdir("../GenGrid_py")
os.system("./gen_grid.py ../../../data/ModelParams/" + site + "/" + modelboxFile + " ../../../data/ModelParams/" + site)
