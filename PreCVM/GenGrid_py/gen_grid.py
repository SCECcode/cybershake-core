#!/usr/bin/env python2

import os
import sys
import math

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

MPI_CMD = config.getProperty("MPI_CMD")
 
def genGridfile(site, outputFile, xlen, ylen, zlen, spacing):
    '''Takes a site, path to output directory, and data to produce gridfile_<site> with the X, Y, and Z lengths and spacing'''
    output = open(outputFile, "w")
    output.write("xlen=%.6f\n" % xlen)
    output.write("    0.0000   %.6f  %e\n" % (xlen, spacing))
    output.write("ylen=%.6f\n" % ylen)
    output.write("    0.0000   %.6f  %e\n" % (ylen, spacing))
    output.write("zlen=%.6f\n" % zlen)
    output.write("    0.0000    %.6f  %e\n" % (zlen, spacing))
    output.flush()
    output.close()

def genBoundfile(gridout, coordfile, boundfile):
    '''Takes a path to a gridout_<site> and model_coords_GC_<site> to produce a model_bounds_GC_<site> file, containing a list of all the (x,y) border coordinates.'''
    input = open(gridout)
    gridoutContents = [line.strip() for line in input.readlines()]
    input.close()
    nx = (gridoutContents[1].split("="))[1]
    intNX = int(nx)
    ny = (gridoutContents[1+intNX+2].split("="))[1]
    intNY = int(ny)
    input = open(coordfile)
    output = open(boundfile, "w")
    for line in input:
        pieces = line.split()
        if int(pieces[2])==0 or int(pieces[2])==intNX-1 or int(pieces[3])==0 or int(pieces[3])==intNY-1:
            output.write(line)
    
    input.close()
    output.flush()
    output.close()


def genGrid(modelboxFile, gridfile, gridout, coordfile, paramfile, boundsfile, freq, sp, gpu=False, depth=-1.0):
	'''Replaces the gen_grid.csh script;  produces a regular grid from a modelbox file.'''
	SPACING = 0.1/freq
	if sp>0:
		SPACING = sp

	#Changed to ZLEN = 50.4 for central CA, since we're propagating over a larger distance, requires GPU or CPU counts get all off
	if gpu:
		ZLEN = 50.4
		#ZLEN = 40.0
	else:
		#ZLEN = 50.0
		ZLEN = 40.0

	if depth>0:
		ZLEN = depth
	#Make sure it is evenly divisible, and nz is even
	#Use ceil because we want the depth to be at least ZLEN
	nz = int(math.ceil(ZLEN/SPACING))
        print("Initial nz is %d\n" % nz)
	if nz*SPACING!=ZLEN:
		ZLEN = nz*SPACING
	if nz%4!=0:
		nz += 4 - (nz%4)
		ZLEN = nz*SPACING
		print("ZLEN adjusted to %f." % ZLEN)
	#if gpu:
	#	#Change ZLEN to 51.2 so it's a multiple of 256 pts
	#	ZLEN = 51.2
	#SPACING = .2

	SPACING = 0.1/freq
	if sp>0:
		SPACING = sp

	print "SPACING = %f, freq = %f, sp= %f\n" % (SPACING, freq, sp)
	modelboxInput = open(modelboxFile)
	modelboxData = [line.strip() for line in modelboxInput.readlines()]
	modelboxInput.close()
    
	site = modelboxData[0]
	# skipping centroid info and parameters header
	modelParamsLine = modelboxData[4]
	# this line is
	#    mlon= <lon>  mlat= <lat>  mrot= <rot>  xlen = <xdist>  ylen= <ydist>
	modelParamsData = modelParamsLine.split()
	model_lon = float(modelParamsData[1])
	model_lat = float(modelParamsData[3])
	model_rot = float(modelParamsData[5])
	xlen = float(modelParamsData[7])
	ylen = float(modelParamsData[9])

	genGridfile(site, gridfile, xlen, ylen, ZLEN, SPACING)
    
	if MPI_CMD=="aprun":
		executable = "aprun -n 1 %s/PreCVM/GenGrid_py/bin/gen_model_cords" % (config.getProperty("CS_PATH"))
	elif MPI_CMD=="jsrun":
		executable = "jsrun -n 1 %s/PreCVM/GenGrid_py/bin/gen_model_cords" % (config.getProperty("CS_PATH"))    
	elif MPI_CMD=='ibrun':
		executable = "ibrun %s/PreCVM/GenGrid_py/bin/gen_model_cords" % (config.getProperty("CS_PATH"))
	parameters = "geoproj=1 gridfile=%s gridout=%s center_origin=1 do_coords=1 nzout=1 name=%s gzip=0 latfirst=0 modellon=%f modellat=%f modelrot=%f" % (gridfile, gridout, coordfile, model_lon, model_lat, model_rot)
	pipe = "> " + paramfile
	command = executable + " " + parameters + " " + pipe
	exitcode = os.system(command)
	print command
	if exitcode!=0:
		print "Exit with code %d\n" % ((exitcode >> 8) & 0xFF)
		sys.exit((exitcode >> 8) & 0xFF)
	genBoundfile(gridout, coordfile, boundsfile)


def main():
	if len(sys.argv) < 9:
		print "Syntax: gen_grid.py <modelboxFile> <gridfile> <gridout> <coordfile> <paramfile> <boundsfile> <frequency> <spacing> [gpu] [depth]"
		print "Example: gen_grid.py USC.modelbox ModelParams/USC/gridfile_USC ModelParams/USC/gridout_USC ModelParams/USC/model_coords_GC_USC ModelParams/USC/model_params_GC_USC ModelParams/USC/model_bounds_GC_USC 0.5"
		sys.exit(1)
	
	modelboxFile = sys.argv[1]
	gridfile = sys.argv[2]
	gridout = sys.argv[3]
	coordfile = sys.argv[4]
	paramfile = sys.argv[5]
	boundsfile = sys.argv[6]
	frequency = float(sys.argv[7])
	spacing = float(sys.argv[8])
	gpu = False
	depth = -1.0
	index = 9
	while index<len(sys.argv):
		if sys.argv[index]=="gpu":
			gpu = True
		else:
			#It's a depth
			depth = float(sys.argv[index])
		index += 1
	
	genGrid(modelboxFile, gridfile, gridout, coordfile, paramfile, boundsfile, frequency, spacing, gpu=gpu, depth=depth)

if __name__=="__main__":
    main()

