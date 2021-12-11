#!/usr/bin/env python

import os
import sys
import difflib

if (len(sys.argv) != 4):
    print "usage: %s rupdir1 rupdir2 startvar" % (sys.argv[0])
    sys.exit(1)
    
rupdir1= sys.argv[1]
rupdir2 = sys.argv[2]
startvar = int(sys.argv[3])

ip = open("%s/%s" % (rupdir1, "variations.list"), "r")
vars = ip.readlines()
ip.close()
count=0
print "Found %d variations" % (len(vars))
for i in xrange(startvar, len(vars)):
    r = vars[i]
    src,rup,numslip,numhypo=r.split()
    for i in xrange(0, int(numslip)):
        for j in xrange(0, int(numhypo)):
            infile="%s/%d/%d/%d/%d_%d.txt.variation-s%04d-h%04d" % \
                    (rupdir1, int(src), int(rup), i, int(src), \
                     int(rup), i, j)
            outfile="%s/%d/%d/%d/%d_%d.txt.variation-s%04d-h%04d" % \
                    (rupdir2, int(src), int(rup), i, int(src), \
                     int(rup), i, j)
            #print "Comparing %s and %s" % (infile, outfile)
            ip1 = open(infile, "r")
            ip2 = open(outfile, "r")
            s = difflib.SequenceMatcher(None, ip1.read(), ip2.read())
            ip1.close()
            ip2.close()
            if (s.ratio() < 1.0):
                print "Mismatch: %s and %s" % (infile, outfile)
            #os.system("/usr/bin/diff -q %s %s" % (infile, outfile))
            if ((count % 1000 == 0) and (count > 0)):
                print "Processed %d variations" % (count)
            count = count + 1;
            #sys.exit(0)

print "Done."
sys.exit(0)
