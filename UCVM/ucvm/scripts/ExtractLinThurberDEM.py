#!/usr/bin/env python
##########################################################
#
# Script: ExtractLinThurberDEM.py
#
# Description: Parse USGS DEMs in GridFloat format to
# produce Lin-Thurber 1km DEM. Topo map is x-y grid
# of 4-byte floats representing elevation in meters.
#
##########################################################


# Basic modules
import os
import sys
import struct
import pyproj
import numpy
from ParseConfig import *
from Interpolate import *
from DEM import *

# DEM Constants
DEM_DIR = '../bin'
DEM_SCRIPT = 'run_dem.sh'


class ExtractLinThurberDEM:
    def __init__(self, conf, neddir, bathdir, outfile):
        self.valid = False
        self.conf = conf
        self.neddir = neddir
        self.bathdir = bathdir
        self.outfile = outfile
        self.valid = True


    def isValid(self):
        return self.valid


    def cleanup(self):
        return


    def _parseModel(self):
        fp = open(self.conf, 'r')
        data = fp.readlines()
        fp.close()

        p = ParseConfig(data)
        p.showDict()
        config = p.getDict()
 
        self.xi = config['proj_xi'].split(',')
        self.yi = config['proj_yi'].split(',')
        for i in xrange(0, 4):
            self.xi[i] = float(self.xi[i])
            self.yi[i] = float(self.yi[i])

        self.size = config['proj_size'].split(',')
        for i in xrange(0, 2):
            self.size[i] = float(self.size[i])

        self.num_z = int(config['num_z'])
        self.spacing = float(config['spacing_dem'])
        self.dims = [int(self.size[0]/self.spacing) + 1, \
                         int(self.size[1]/self.spacing) + 1]
        self.z_vals = config['z_vals'].split(',')
        for i in xrange(0, self.num_z):
            self.z_vals[i] = float(self.z_vals[i])

        return(0)


    def _getGrid(self):
        grid = []
        i = Interpolate(self.xi, self.yi, self.size)
        for y in xrange(0, self.dims[1]):
            for x in xrange(0, self.dims[0]):
                lon,lat = i.binterp(x * self.spacing, \
                                        y * self.spacing, inverse=True)
                grid.append([lon, lat])

        return(grid);


    def _queryDEM(self, points):
        infile = 'dem.in'
        outfile = 'dem.out'
        fp = open(infile, 'w')
        for p in points:
            fp.write("%lf %lf\n" % (p[0], p[1]))
        fp.close()

        dem = DEM(DEM_DIR, DEM_SCRIPT, self.neddir, self.bathdir)
        dem.run(infile, outfile)

        fp = open(outfile, 'r')
        demdata = fp.readlines()
        fp.close()

        os.remove(infile)
        os.remove(outfile)

        data = numpy.arange(self.dims[0]*self.dims[1], dtype=float).reshape(self.dims[1], self.dims[0])
        for y in xrange(0, self.dims[1]):
            for x in xrange(0, self.dims[0]):
                tokens = demdata[y*self.dims[0]+x].split()
                elev = float(tokens[3])
                flag = int(tokens[4])
                if (flag):
                    data[y][x] = elev
                else:
                    print "No elev data for point %lf, %lf" % \
                        (points[y*self.dims[0]+x][0], \
                             points[y*self.dims[0]+x][1])
                    data[y][x] = 0.0

        return(data)


    def main(self):

        # Parse LT model
        print "Parsing model config"
        self._parseModel()

        # Generate grid with bilinear interp
        print "Generating 2D grid"
        points = self._getGrid()

        # Query DEMs
        print "Querying DEM"
        data = self._queryDEM(points)

        # Write topo file, x-fast
        print "Writing DEM %s" % (self.outfile)
        fp = open(self.outfile, 'wb')
        for y in xrange(0, self.dims[1]):
            for x in xrange(0, self.dims[0]):
                val = struct.pack('f', data[y][x])
                fp.write(val)
        fp.close()
        return 0


def usage():
    print "usage: %s <lt_conf> <neddir> <bathdir> <outfile>" % (sys.argv[0])
    return


if __name__ == '__main__':
    if (len(sys.argv) != 5):
        usage()
        sys.exit(1)

    conf = sys.argv[1]
    neddir = sys.argv[2]
    bathdir = sys.argv[3]
    outfile = sys.argv[4]

    prog = ExtractLinThurberDEM(conf, neddir, bathdir, outfile)
    sys.exit(prog.main())
