#!/usr/bin/env python

import sys
import os

if len(sys.argv)<3:
	print "Usage: %s <1D velocity model> <output>" % sys.argv[0]
	sys.exit(1)

velocity_model = sys.argv[1]
out_filename = sys.argv[2]
c0=57
c1=34
i=0
vm_fp = open(velocity_model, "r")
vel_in_data = vm_fp.readlines()
vm_fp.close()
#vel_out = "local_%s" % velocity_model
#vel_out_fp = open(vel_out, "w")
vel_out_fp = open(out_filename, "w")
for line in vel_in_data:
	i += 1
        if line.startswith("#") or line.startswith("%"):
                continue
        pieces = line.split()
        if len(pieces)>=4:
                th=float(pieces[0])
                vp=float(pieces[1])
                vs=float(pieces[2])
                dn=float(pieces[3])
                qs=c0+c1*vs
                if i==len(vel_in_data):
                        th=0.0
                vel_out_fp.write("%8.4f %8.4f %8.4f %8.4f %8.2f %8.2f\n" % (th, vp, vs, dn, qs, qs))
        else:
                vel_out_fp.write(line)
vel_out_fp.flush()
vel_out_fp.close()
sys.exit(0)
