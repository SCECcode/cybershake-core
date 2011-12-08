#!/usr/bin/env python

import os
import sys
import glob
import shutil

indir = sys.argv[1]
outdir = sys.argv[2]

print indir, outdir

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
    i = i + 1
    if (i % 500 == 0):
        print "Processed %d files" % (i)

sys.exit(0)
