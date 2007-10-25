#!/usr/bin/python

import os
import sys
 
def genGridfile(site, outputDirectory, xlen, ylen, zlen, spacing):
    '''Takes a site, path to output directory, and data to produce gridfile_<site> with the X, Y, and Z lengths and spacing'''
    gridfileName = outputDirectory + "/gridfile_" + site
    output = open(gridfileName, "w")
    output.write("xlen=%.6f\n" % xlen)
    output.write("    0.0000   %.4f  %e\n" % (xlen, spacing))
    output.write("ylen=%.6f\n" % ylen)
    output.write("    0.0000   %.4f  %e\n" % (ylen, spacing))
    output.write("zlen=%.6f\n" % zlen)
    output.write("    0.0000    %.4f  %e\n" % (zlen, spacing))
    output.flush()
    output.close()
    return gridfileName


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


def genGrid(modelboxFile, outputDirectory):
    '''Replaces the gen_grid.csh script;  produces a regular grid from a modelbox file.'''
    ZLEN = 40.0
    SPACING = .2

    modelboxInput = open(modelboxFile)
    modelboxData = [line.strip() for line in modelboxInput.readlines()]
    modelboxInput.close()

    site = modelboxData[0]
    # skipping centroid info and parameters header
    modelParamsLine = modelboxData[4]
    # this line is
    #    mlon= <lon>  mlat= <lat>  mrot= <rot>  xlen = <xdist>  ylen= <ydistdfgCVS/CVS;
    modelParamsData = modelParamsLine.split()
    model_lon = float(modelParamsData[1])
    model_lat = float(modelParamsData[3])
    model_rot = float(modelParamsData[5])
    xlen = float(modelParamsData[7])
    ylen = float(modelParamsData[9])

    gridfileName = genGridfile(site, outputDirectory, xlen, ylen, ZLEN, SPACING)
    gridoutName = outputDirectory + "/gridout_" + site
    coordfileName = outputDirectory + "/model_coords_GC_" + site
    paramfileName = outputDirectory + "/model_params_GC_" + site
    boundsfileName = outputDirectory + "/model_bounds_GC_" + site

    executable = "bin/gen_model_cords"
    parameters = "geoproj=1 gridfile=%s gridout=%s center_origin=1 do_coords=1 nzout=1 name=%s gzip=0 latfirst=0 modellon=%f modellat=%f modelrot=%f" % (gridfileName, gridoutName, coordfileName, model_lon, model_lat, model_rot)
    pipe = "> " + paramfileName
    command = executable + " " + parameters + " " + pipe
    print command

    os.system(command)
    genBoundfile(gridoutName, coordfileName, boundsfileName)



def main():
    if len(sys.argv) < 3:
        print "Syntax: gen_grid.py <modelboxFile> <outputDirectory>"
        print "Example: gen_grid.py USC.modelbox ModelParams/"
        sys.exit()
	
    modelboxFile = sys.argv[1]
    outputDirectory = sys.argv[2]
    genGrid(modelboxFile, outputDirectory)

if __name__=="__main__":
    main()

