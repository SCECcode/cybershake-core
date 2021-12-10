#!/usr/bin/env python

import sys
import os
import config

def build_media(site, gridout, rwg_vel_prefix, media_out):
	cs_path = config.getProperty("CS_PATH")
	
	#determine NX, NY, NZ
	fp_in = open(gridout, "r")
	data = fp_in.readlines()
	fp_in.close()
	nx = int((data[1].split("="))[1])
	ny = int((data[1+nx+2].split("="))[1])
	nz = int((data[1+nx+2+ny+2].split("="))[1])
	
	cmd = "%s/SgtHead/bin/reformat_velocity %d %d %d %s %s" % (cs_path, nx, ny, nz, rwg_vel_prefix, media_out)
	print(cmd)
	exitcode = os.system(cmd)
	return ((exitcode >> 8) & 0xFF)

