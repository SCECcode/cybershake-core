#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
import os
import sys
from pylab import *
import PlotSRF

bin_path = "/work/00940/tera3d/CyberShake/software/SeisPSA_header_buf/SlipModel/StandRupFormat"

srf_file = sys.argv[1]
slipfile = "%s.slip" % (srf_file.split(".")[0])
tinit_file = "%s.tinit" % (srf_file.split(".")[0])

cmd = "%s/srf2xyz calc_xy=0 type=slip nseg=-1 infile=%s > %s" % (bin_path, srf_file, slipfile)
print cmd
os.system(cmd)
cmd = "%s/srf2xyz calc_xy=0 type=tinit nseg=-1 infile=%s > %s" % (bin_path, srf_file, tinit_file)
print cmd
os.system(cmd)

plotter = PlotSRF.PlotSRF()
plotter.plot(os.path.basename(srf_file), srf_file, "slip", ".")


