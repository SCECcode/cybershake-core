#!/bin/env python
"""
Plots slip distribution for a SRF
$Id: PlotSRF.py 992 2012-09-17 19:04:41Z fsilva $
"""

# Import Python modules
import os
import sys
import numpy as np
import matplotlib
if (matplotlib.get_backend() != 'agg'):
    matplotlib.use('Agg') # Disables use of Tk/X11
import matplotlib.colors as mcolors 
import matplotlib.cm as cm
from pylab import *

# Slip range in cm
SLIP_X_FACTOR = 20.0
SLIP_Y_FACTOR = 5.0

# Contour intervals
CONTOUR_INTERVALS = 10.0

class PlotSRF:
    def __init__(self):
        """
        Initializes the PlotSRF class with None values
        """
        self.srf_file = None
        self.fault_len = None
        self.fault_width = None
        self.dim_len = None
        self.dim_wid = None

    def read_srf(self, file, numx, numy):
        """
        Read in fault file
        """
        ip = open(file)
        slips = ip.readlines()
        ip.close()

        data = np.arange(numx * numy, dtype=float).reshape(numy, numx)

        # Data is x-fast
        for y in xrange(0, numy):
            for x in xrange(0, numx):
		#skip first line
                tokens = slips[y * (numx) + x + 1].split()
                data[y][x] = tokens[2]

        return(data)

    def get_srf_params(self):
        """
        Reads fault_len, width, dlen, and dwid from the srd file
        """
        srf_params = None
        srf = open(self.srf_file, 'r')
        for line in srf:
            if line.startswith("PLANE 1"):
                # Found the plane line, the next one should have what we need
                srf_params = srf.next()
                break
        srf.close()
        if srf_params is None:
            print ("Cannot determine parameters from SRF file %s" %
                   self.srf_file)
            sys.exit(1)
        srf_params = srf_params.strip()
        srf_params = srf_params.split()
        # Make sure we have the correct number of pieces
        if len(srf_params) != 6:
            print ("Cannot parse params from SRF file %s" %
                   self.srf_file)
            sys.exit(1)
        self.dim_len = int(srf_params[2])
        self.dim_wid = int(srf_params[3])
        self.fault_len = float(srf_params[4])
        self.fault_width = float(srf_params[5])

    def plot(self, plottitle, srffile, plottype, outdir):
        """
        Produce the plot
        """
        self.srf_file = srffile
        # Get SRF parameters
        self.get_srf_params()

        # Plot dimensions
        self.dims = [self.dim_len, self.dim_wid]
        extents = [0.0, self.fault_len, self.fault_width, 0.0]

        # Read in SRF slips
        slipfile = "%s.slip" % (os.path.splitext(srffile)[0])
        slips = self.read_srf(slipfile, self.dims[0], self.dims[1])

        # Read in SRF tinits
        tinitfile = "%s.tinit" % (os.path.splitext(srffile)[0])
        tinits = self.read_srf(tinitfile, self.dims[0], self.dims[1])

        # Find avg/max slip
        avgslip = 0.0
        minslip = 100000.0
        maxslip = 0.0
        for y in xrange(0, self.dims[1]):
            for x in xrange(0, self.dims[0]):
                if (slips[y][x] > maxslip):
                    maxslip = slips[y][x]
                if (slips[y][x] < minslip):
                    minslip = slips[y][x]
                avgslip = avgslip + slips[y][x]
        avgslip = avgslip / (self.dims[0] * self.dims[1])

        # Set plot dims
        gcf().set_size_inches(6, 8)
        gcf().clf()

        # Set title and adjust title y-position
        t = title("%s\nAvg/Max Slip = %d/%d" % (plottitle,
                                                int(avgslip),
                                                int(maxslip)), size=12)
        t.set_y(1.05) 

        # Setup slip color scale
        cmap = cm.hot_r
        d = int(maxslip / SLIP_X_FACTOR + 0.0)
        while (SLIP_X_FACTOR * d < 0.9 * maxslip):
            d = d + 1
        colormin = 0.0
        colormax = float(SLIP_X_FACTOR * d)
        colorint = float(SLIP_Y_FACTOR * d)
        norm = mcolors.Normalize(vmin=colormin, vmax=colormax)

        # Plot slips
        imshow(slips, cmap=cmap, norm=norm, extent=extents,
               interpolation='nearest')

        # Freeze the axis extents
        gca().set_autoscale_on(False)
        xlabel("Along Strike (km)", size=8)
        ylabel("Down Dip (km)", size=8)

        # Set font size
        for tick in gca().get_xticklabels():
            tick.set_fontsize(8)
        for tick in gca().get_yticklabels():
            tick.set_fontsize(8)

        # Setup slip color scale
        cb = colorbar(orientation='horizontal', shrink=0.5,
                      ticks=linspace(colormin, colormax,
                                     (colormax/colorint) + 1))
        cb.set_label('Slip (cm)', fontsize=8)
        for tick in cb.ax.get_xticklabels():
            tick.set_fontsize(8)

        # Setup tinit contours
        mintinit = 100000.0
        maxtinit = 0.0
        for y in xrange(0, self.dims[1]):
            for x in xrange(0, self.dims[0]):
                if (tinits[y][x] > maxtinit):
                    maxtinit = tinits[y][x]
                if (tinits[y][x] < mintinit):
                    mintinit = tinits[y][x]
        tinit_interval = (maxtinit - mintinit) / 10.0
 
        # Plot tinit contours
        contour(tinits,
                linspace(mintinit, maxtinit, CONTOUR_INTERVALS + 1.0),
                origin='upper', extent=extents, colors = 'k')

        outfile = os.path.join(outdir,
                               "%s.png" %
                               (os.path.splitext(srffile)[0]))
        print "Saving plot to %s" % (outfile)
        savefig(outfile, format="png", transparent=False)
        return

def usage():
    print ("usage: %s <srf file> <slip|velocity> <outdir>" %
           (sys.argv[0]))
    return

if __name__ == '__main__':
    if (len(sys.argv) != 4):
        usage()
        sys.exit(1)

    srffile = sys.argv[1]
    plottype = sys.argv[2]
    outdir = sys.argv[3]

    plottitle = 'Rupture Model for %s' % (os.path.basename(srffile))

    plotter = PlotSRF()
    plotter.plot(plottitle, srffile, plottype, outdir)
    sys.exit(0)
   
