#!/usr/bin/env python

import sys
import os
import optparse

parser = optparse.OptionParser()
parser.add_option("--site", dest="site", action="store", help="Site name")
parser.add_option("--erf_id", dest="erf_id", action="store", type="int", help="ERF ID")
parser.add_option("--modelbox", dest="modelbox", action="store", help="Path to modelbox file (output)")
parser.add_option("--gridfile", dest="gridfile", action="store", help="Path to gridfile (output)")
parser.add_option("--gridout", dest="gridout", action="store", help="Path to gridout (output)")
parser.add_option("--coordfile", dest="coordsfile", action="store", help="Path to coorfile (output)")
parser.add_option("--paramsfile", dest="paramsfile", action="store", help="Path to paramsfile (output)")
parser.add_option("--boundsfile", dest="boundsfile", action="store", help="Path to boundsfile (output)")
parser.add_option("--frequency", dest="frequency", action="store", type="float", help="Frequency")
parser.add_option("--gpu", dest="gpu_arg", action="store_true", default=False, help="Use GPU box settings.")
parser.add_option("--spacing", dest="spacing", action="store", type="float", help="Override default spacing with this value.")
parser.add_option("--server", dest="server", action="store", default="focal.usc.edu", help="Address of server to query in creating modelbox, default is focal.usc.edu.")

(option, args) = parser.parse_args()

site = option.site
erfID = option.erf_id
modelbox = os.path.abspath(option.modelbox)
gridfile = os.path.abspath(option.gridfile)
gridout = os.path.abspath(option.gridout)
coordsfile = os.path.abspath(option.coordsfile)
paramsfile = os.path.abspath(option.paramsfile)
boundsfile = os.path.abspath(option.boundsfile)

if site==None or erfID==None or modelbox==None or gridfile==None or gridout==None or coordsfile==None or paramsfile==None or boundsfile==None:
	print "All of site, ERF ID, modelbox, gridfile, gridout, coordsfile, paramsfile, and boundsfile must be specified."
	parser.print_help()
	sys.exit(-1)

frequency = 0.5
if option.frequency is not None:
	frequency = option.frequency
use_gpu = option.gpu_arg
spacing = -1.0
if option.spacing is not None:
	spacing = option.spacing
server = option.server

os.chdir(os.path.join(sys.path[0], "Modelbox"))
gpu_arg = ""
if use_gpu:
	gpu_arg = "gpu"
 
exitcode = os.system("./get_modelbox.py %s %s %s %f %s %s" % (site, erfID, modelbox, spacing, server, gpu_arg))
if exitcode!=0:
	sys.exit((exitcode >> 8) & 0xFF)
os.chdir("../GenGrid_py")
exitcode = os.system("./gen_grid.py %s %s %s %s %s %s %f %f %s" % (modelbox, gridfile, gridout, coordsfile, paramsfile, boundsfile, frequency, spacing, gpu_arg))
sys.exit((exitcode >> 8) & 0xFF)
