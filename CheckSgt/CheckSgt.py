#!/usr/bin/env python

import sys
import os
import md5


# Global vars
sgtfile = ""
md5file = ""


def init():
    global sgtfile
    global md5file
        
    # Get number of command-line arguments
    argc = len(sys.argv)
    
    # Parse command line arguments
    if (argc < 3):
        print "Usage: " + sys.argv[0] + " <sgt file> <md5 file>"
        print "Example: " + sys.argv[0] + " USC_fx.sgt USC_fx.sgt.md5"
        return 1
            
    sgtfile = sys.argv[1]
    md5file = sys.argv[2]

    print "Configuration:"
    print "SGT File:\t" + sgtfile
    print "MD5 File:\t" + md5file + "\n"
    
    # Check that the files exist
    if (not os.path.isfile(sgtfile)):
        print "SGT file " + sgtfile + " not found"
        return 1
    if (not os.path.isfile(md5file)):
        print "MD5 file " + md5file + " not found"
        return 1  
    
    return 0


def main():
    # Load the sum from the saved md5 file
    oldmd5val = ""
    try:
        oldmd5file = open(md5file, 'r')
        line = oldmd5file.readline()
        oldmd5val = line.split(" ")[0]
        oldmd5file.close()
        if (len(oldmd5val) != 32):
            print "Invalid md5sum found: " + oldmd5val
            return 1
    except:
        print "Unable to read " + md5file
        return 1
    
    print "Old md5sum: " + oldmd5val
    
    # Read in the SGT file and compute new md5
    m = md5.new()
    sgt = open(sgtfile, 'r') # open in binary mode
    while True:
        buf = sgt.read(1024)
        if len(buf) == 0:
            break # end of file
        m.update(buf)

    newmd5val = m.hexdigest()
    print "New md5sum: " +  newmd5val

    # Compare the old and new md5 values
    if (oldmd5val != newmd5val):
        print "md5 checksums do not match!"
        return 1
    
    return 0


def cleanup():
    return 0


if __name__ == '__main__':
    if (init() != 0):
        sys.exit(1)
    if (main() != 0):
        sys.exit(1)
    cleanup()
    sys.exit(0)
