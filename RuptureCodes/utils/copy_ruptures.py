#!/usr/bin/env python

import os
import sys
import glob
import shutil


ip = open("%s/%s" % (sys.argv[1], "ruptures.list"), "r")
ruptures = ip.readlines()
ip.close()
for r in ruptures:
    src,rup,mag,file=r.split()
    rupdir="%s/%s/%s" % (sys.argv[2], src, rup)
    try:
        os.makedirs(rupdir)
    except:
        pass
    infile="%s%s" % (sys.argv[1], file)
    outfile="%s%s" % (sys.argv[2], file)
    print "Copying %s to %s" % (infile, outfile)
    shutil.copy2(infile, outfile)


shutil.copy2("%s/%s" % (sys.argv[1], "ruptures.list"), \
             "%s/%s" % (sys.argv[2], "ruptures.list"))
