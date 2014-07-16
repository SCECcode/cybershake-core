#!/usr/bin/env python

'''This script constructs an IN3D file for a CyberShake AWP-ODC-SGT run.
Some values are taken from data/IN3D.ref;  the others are determined based on the run.
X and Y are swapped from the RWG system, which is what all the input files are based on.
Since this is the GPU version, we have to make sure NX/PX and NY/PY are both even.'''

import sys
import os
import math

def build_IN3D(site, gridout, awp_comp, frequency, proc):
	fp_in = open("%s/data/IN3D.ref" % (sys.path[0]), "r")
	data = fp_in.readlines()
	fp_in.close()

	param = dict()
	
	for line in data:
		pieces = line.split("=")
		param[pieces[0]] = pieces[1].strip()

	#determine igreen
	igreen = 0
	if awp_comp=='x':
		#swap x and y
		comp = 'y'
		igreen = 4
	elif awp_comp=='y':
		#swap x and y
		comp = 'x'
		igreen = 5
	elif awp_comp=='z':
		comp = 'z'
		igreen = 6
	else:
		print "Component %s not recognized, aborting." % comp
		sys.exit(2)
	
	param["igreen"] = igreen

	#determine DH, DT, NST, READ_STEP, WRITE_STEP, FP
	param["DH"] = round(100.0/frequency, 1)
	param["DT"] = 0.005/frequency
	param["NST"] = int(frequency*40000.0)
	param["FP"] = frequency
	param["READ_STEP"] = param["NST"]
	#Divide by 10 because WRITE_STEP is in units of # of steps being written, not total # of timesteps
	#So if 20000 simulated timesteps and decimation of 10, WRITE_STEP = 2000
	param["WRITE_STEP"] = param["READ_STEP"]/int(param["NTISKP_SGT"])
	param["WRITE_STEP2"] = param["WRITE_STEP"]
	
	#determine NX, NY, NZ
	fp_in = open(gridout, "r")
	data = fp_in.readlines()
	fp_in.close()
	#remember, X and Y are flipped
	ny = int((data[1].split("="))[1])
	nx = int((data[1+ny+2].split("="))[1])
	nz = int((data[1+ny+2+nx+2].split("="))[1])
	
	param["NX"] = nx
	param["NY"] = ny
	param["NZ"] = nz

	#Check proc values to make sure they are evenly divisible
	if nx % proc[0] != 0:
		print "PX %d must be a factor of NX %d, aborting." % (proc[0], nx)
		sys.exit(2)
        if ny % proc[1] != 0:
                print "PY %d must be a factor of NY %d, aborting." % (proc[1], ny)
                sys.exit(2)
        if nz % proc[2] != 0:
                print "PZ %d must be a factor of NZ %d, aborting." % (proc[2], nz)
                sys.exit(2)	
	
	param["NPX"] = proc[0]
	param["NPY"] = proc[1]
	param["NPZ"] = proc[2]
	
	#doesn't really matter, but set N<BG|ED><comp> values
	param["NBGX"] = 1
	param["NEDX"] = nx
	param["NBGY"] = 1
	param["NEDY"] = ny
	param["NBGZ"] = 1
	param["NEDZ"] = nz
	param["NBGX2"] = 1
	param["NEDX2"] = nx
	param["NBGY2"] = 1
	param["NEDY2"] = ny
	param["NBGZ2"] = 1
	param["NEDZ2"] = nz
	
	#paths to INSRC, INVEL, INSGT, SGTGR0
	param["INSRC"] = "comp_%s/input/%s_f%s_src" % (awp_comp, site, awp_comp)
	param["INVEL"] = "comp_%s/input/awp.%s.media" % (awp_comp, site)
	param["INSGT"] = "comp_%s/input/awp.%s.cordfile" % (awp_comp, site)
	param["SGTGRO"] = "comp_%s/output_sgt/awp-strain-%s-f%s" % (awp_comp, site, awp_comp)

	
	#write IN3D file
	fp_out = open("IN3D.%s.%s" % (site, awp_comp), "w")
	fp_out.write("%9s igreen\n" % (param["igreen"]))
	fp_out.write("%9s TMAX\n\n" % (param["TMAX"]))
	fp_out.write("%9s DH\n" % (param["DH"]))
	fp_out.write("%9s DT\n\n" % (param["DT"]))
	fp_out.write("%9s NPC\n\n" % (param["NPC"]))
	fp_out.write("%9s ND\n" % (param["ND"]))
	fp_out.write("%9s ARBC\n" % (param["ARBC"]))
	fp_out.write("%9s PHT\n\n" % (param["PHT"]))
	fp_out.write("%9s NSRC\n" % (param["NSRC"]))
	fp_out.write("%9s NST\n\n" % (param["NST"]))
	fp_out.write("%9d NX\n" % (param["NX"]))
	fp_out.write("%9d NY\n" % (param["NY"]))
	fp_out.write("%9d NZ\n\n" % (param["NZ"]))
	fp_out.write("%9d NPX\n" % (param["NPX"]))
	fp_out.write("%9d NPY\n" % (param["NPY"]))
	fp_out.write("%9d NPZ\n\n" % (param["NPZ"]))
	fp_out.write("%9s IFAULT\n" % (param["IFAULT"]))
	fp_out.write("%9s CHECKPOINT\n" % (param["CHECKPOINT"]))
	fp_out.write("%9s ISFCVLM\n" % (param["ISFCVLM"]))
	fp_out.write("%9s IMD5\n" % (param["IMD5"]))
	fp_out.write("%9s IVELOCITY\n" % (param["IVELOCITY"]))
	fp_out.write("%9s MEDIARESTART\n" % (param["MEDIARESTART"]))
	fp_out.write("%9s NVAR\n" % (param["NVAR"]))
	fp_out.write("%9s IOST\n" % (param["IOST"]))
	fp_out.write("%9s PARTDEG\n" % (param["PARTDEG"]))
	fp_out.write("%9s IO_OPT\n" % (param["IO_OPT"]))
	fp_out.write("%9s PERF_MEAS\n" % (param["PERF_MEAS"]))
	fp_out.write("%9s IDYNA\n" % (param["IDYNA"]))
	fp_out.write("%9s SOCALQ\n\n" % (param["SOCALQ"]))
	fp_out.write("%9s NVE\n\n" % (param["NVE"]))
	fp_out.write("%9s MU_S\n" % (param["MU_S"]))
	fp_out.write("%9s MU_D\n\n" % (param["MU_D"]))
	fp_out.write("%9s FL\n" % (param["FL"]))
	fp_out.write("%9s FH\n" % (param["FH"]))
	fp_out.write("%9s FP\n\n" % (param["FP"]))
	fp_out.write("%9s READ_STEP\n" % (param["READ_STEP"]))
	fp_out.write("%9s WRITE_STEP\n" % (param["WRITE_STEP"]))
	fp_out.write("%9s WRITE_STEP2\n\n" % (param["WRITE_STEP2"]))
	fp_out.write("%9s NBGX\n" % (param["NBGX"]))
	fp_out.write("%9s NEDX\n" % (param["NEDX"]))
	fp_out.write("%9s NSKPX\n" % (param["NSKPX"]))
	fp_out.write("%9s NBGY\n" % (param["NBGY"]))
	fp_out.write("%9s NEDY\n" % (param["NEDY"]))
	fp_out.write("%9s NSKPY\n" % (param["NSKPY"]))
	fp_out.write("%9s NBGZ\n" % (param["NBGZ"]))
	fp_out.write("%9s NEDZ\n" % (param["NEDZ"]))
	fp_out.write("%9s NSKPZ\n\n" % (param["NSKPZ"]))
	#NTISKP needs to be larger than # of timesteps so that no velocity data is written
	fp_out.write("%9s NTISKP\n" % (2*int(param["NST"])))
	fp_out.write("%9s NBGX2\n" % (param["NBGX2"]))
	fp_out.write("%9s NEDX2\n" % (param["NEDX2"]))
	fp_out.write("%9s NSKPX2\n" % (param["NSKPX2"]))
	fp_out.write("%9s NBGY2\n" % (param["NBGY2"]))
	fp_out.write("%9s NEDY2\n" % (param["NEDY2"]))
	fp_out.write("%9s NSKPY2\n" % (param["NSKPY2"]))
	fp_out.write("%9s NBGZ2\n" % (param["NBGZ2"]))
	fp_out.write("%9s NEDZ2\n" % (param["NEDZ2"]))
	fp_out.write("%9s NSKPZ2\n\n" % (param["NSKPZ2"]))
	fp_out.write("%9s NTISKP2\n\n" % (param["NTISKP"]))
        fp_out.write("%9s NTISKP_SGT\n\n" % (param["NTISKP_SGT"]))
	fp_out.write("'comp_%s/%s' CHKP\n" % (awp_comp, param["CHKP"]))
	fp_out.write("'comp_%s/%s' CHKJ\n\n" % (awp_comp, param["CHKJ"]))
	fp_out.write("'%s' INSRC\n" % (param["INSRC"]))
	fp_out.write("'%s' INVEL\n\n" % (param["INVEL"]))
	fp_out.write("'%s' INSGT\n\n" % (param["INSGT"]))
	fp_out.write("'comp_%s/%s' SXRGO\n" % (awp_comp, param["SXRGO"]))
	fp_out.write("'comp_%s/%s' SYRGO\n" % (awp_comp, param["SYRGO"]))
	fp_out.write("'comp_%s/%s' SZRGO\n\n" % (awp_comp, param["SZRGO"]))
	fp_out.write("'comp_%s/%s' SXRGO2\n" % (awp_comp, param["SXRGO2"]))
	fp_out.write("'comp_%s/%s' SYRGO2\n" % (awp_comp, param["SYRGO2"]))
	fp_out.write("'comp_%s/%s' SZRGO2\n\n" % (awp_comp, param["SZRGO2"]))
	fp_out.write("'%s' SGTGRO\n\n" % (param["SGTGRO"]))
	fp_out.write("'comp_%s/%s' SGSN\n" % (awp_comp, param["SGSN"]))
	
	fp_out.flush()
	fp_out.close()

	return 0
