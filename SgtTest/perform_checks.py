#!/usr/bin/env python

'''This performs a check of SGT size, then calls the nan check code.'''

import sys
import os
import struct

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

if len(sys.argv)<3:
	print "Usage: %s <SGT file> <SGT header file>" % sys.argv[0]
	sys.exit(1)

sgt_filename = sys.argv[1]
sgt_header = sys.argv[2]

#Get actual file size
sgt_filesize = os.path.getsize(sgt_filename)

#Read out header information in header file to check for match
'''
struct sgtmaster
   {
   int geoproj;     /* =0: RWG local flat earth; =1: RWG great circle arcs; =2: UTM *
/
   float modellon;  /* longitude of geographic origin */
   float modellat;  /* latitude of geographic origin */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float xshift;    /* xshift of cartesian origin from geographic origin */
   float yshift;    /* yshift of cartesian origin from geographic origin */
   int globnp;      /* total number of SGT locations (entire model) */
   int localnp;     /* local number of SGT locations (this file only) */
   int nt;          /* number of time points                                */
   };
'''
with open(sgt_header) as fp_in:
	sgtmaster_str = fp_in.read(36)
	globnp = struct.unpack('i', sgtmaster_str[24:28])[0]
	nt = struct.unpack('i', sgtmaster_str[32:36])[0]
	fp_in.close()


FLOAT_SIZE = 4
SGT_COMPONENTS = 6

expected_size = globnp*nt*SGT_COMPONENTS*FLOAT_SIZE

if expected_size!=sgt_filesize:
	print "Error: sgtmaster header file %s leads us to expect %d points and nt=%d, for an SGT filesize of %d, but the SGT file actually has size %d.  Aborting." % (sgt_header, globnp, nt, expected_size, sgt_filesize)
	sys.exit(2)

cs_path = config.getProperty("CS_PATH")
mpi_cmd = config.getProperty("MPI_CMD")

#Make sure it gets executed on a compute node
if mpi_cmd=="aprun":
	mpi_cmd = "aprun -n 1"
elif mpi_cmd=="mpirun":
	mpi_cmd = "mpirun -np 1"

cmd = "%s %s/SgtHead/bin/check_for_nans %s" % (mpi_cmd, cs_path, sgt_filename)
exitcode = os.system(cmd)
if exitcode!=0:
	print "Error when checking for NaNs, zeros."
sys.exit(exitcode)
