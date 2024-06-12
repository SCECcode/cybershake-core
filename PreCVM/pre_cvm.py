#!/usr/bin/env python

import sys
import os
import optparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config


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
parser.add_option("--spacing", dest="spacing", action="store", type="float", help="Override default spacing with this value, in km.")
parser.add_option("--server", dest="server", action="store", default="focal.usc.edu", help="Address of server to query in creating modelbox, default is focal.usc.edu. Can specify sqlite database with sqlite://<file>.  ***Note that if LOCAL_DB is defined in cybershake.cfg file, that value will be used instead of what is given here.***")
parser.add_option("--bounding-box", dest="bbox", action="store_true", default=False, help="Assume (StartLat, StartLon) and (EndLat, EndLon) represent 2 corners of a box, all 4 corners of which must be inside the volume (as opposed to only requiring those 2 points)")
parser.add_option("--tight-box", dest="tight", action="store_true", default=False, help="Use a box with 20 km padding (the default is 50 km)")
parser.add_option("--depth", dest="depth", action="store", type="float", help="Override default depth with this value, in km")
parser.add_option("--rotation", dest="rotation", action="store", type="float", help="Override default rotation with this value")
parser.add_option("--padding", dest="padding", action="store", type="float", help="Override default padding of 50 km with this value, in km")

(option, args) = parser.parse_args()

DEFAULT_ROTATION = -55.0

site = option.site
erfID = option.erf_id
modelbox = os.path.abspath(option.modelbox)
gridfile = os.path.abspath(option.gridfile)
gridout = os.path.abspath(option.gridout)
coordsfile = os.path.abspath(option.coordsfile)
paramsfile = os.path.abspath(option.paramsfile)
boundsfile = os.path.abspath(option.boundsfile)

if site==None or erfID==None or modelbox==None or gridfile==None or gridout==None or coordsfile==None or paramsfile==None or boundsfile==None:
	print("All of site, ERF ID, modelbox, gridfile, gridout, coordsfile, paramsfile, and boundsfile must be specified.")
	parser.print_help()
	sys.exit(-1)

frequency = 0.5
if option.frequency is not None:
	frequency = option.frequency
use_gpu = option.gpu_arg
spacing = 0.1/frequency
if option.spacing is not None:
	spacing = option.spacing
server = option.server
bbox_arg = ""
if option.bbox:
	bbox_arg = "bbox"
tight_arg = ""
if option.tight:
	tight_arg = "tight"
depth = -1.0
if option.depth:
	depth = option.depth
rotation = DEFAULT_ROTATION
if option.rotation is not None:
	rotation = option.rotation
padding = -1.0
if option.padding is not None:
        padding = option.padding

os.chdir(os.path.join(sys.path[0], "Modelbox"))
gpu_arg = ""
if use_gpu:
	gpu_arg = "gpu"

try:
	local_db = config.getProperty('LOCAL_DB')
	server = local_db
except KeyError:
	pass

cmd = "./get_modelbox.py %s %s %s %f %s %f %s %s %s %s" % (site, erfID, modelbox, spacing, server, rotation, padding, gpu_arg, bbox_arg, tight_arg)
print(cmd)
exitcode = os.system(cmd)
if exitcode!=0:
	sys.exit((exitcode >> 8) & 0xFF)

os.chdir("../GenGrid_py")
cmd = "./gen_grid.py %s %s %s %s %s %s %f %f %s %f" % (modelbox, gridfile, gridout, coordsfile, paramsfile, boundsfile, frequency, spacing, gpu_arg, depth)
print(cmd)
exitcode = os.system(cmd)
sys.exit((exitcode >> 8) & 0xFF)
