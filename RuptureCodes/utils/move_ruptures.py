#!/usr/bin/env python

import os
import sys
import glob
import shutil

if (len(sys.argv) != 3):
    print "usage: %s <indir> <outdir>" % (sys.argv[0])
    sys.exit(1)
    
indir = sys.argv[1]
outdir = sys.argv[2]

print indir, outdir

op = open("%s/%s" % (outdir, "ruptures.list"), "w")
if (op == NULL):
    print "Failed to open rupture index file"
    sys.exit(0)
    
ruplist = glob.glob("%s/*.txt" % (indir))
i = 0
for file in ruplist:
    f = os.path.basename(file)
    n, ext = f.split(".")
    src,rup = n.split("_")
    #print file, src, rup
    outpath = "%s/%s/%s" % (outdir, src, rup)
    if (not os.path.exists(outpath)):
        os.makedirs(outpath)
    outfile = "%s/%s" % (outpath, f)
    shutil.move(file, outfile)
    op.write("%s %s %f /%s/%s/%s\n" % (src, rup, 1.0, src, rup, f))
    i = i + 1
    if (i % 500 == 0):
        print "Processed %d files" % (i)

op.close()
sys.exit(0)
