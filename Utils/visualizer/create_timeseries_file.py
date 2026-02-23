#!/usr/bin/env python

import sys
import os
import math

if len(sys.argv)<10:
	print "Usage: %s <X dim in output> <Y dim in output> <grid decimation> <NT in sim> <DT of sim> <NTISKP> <timesteps output per file> <prefix> <output file>" % sys.argv[0]
	sys.exit(1)

x_dim = int(sys.argv[1])
y_dim = int(sys.argv[2])
grid_dec = int(sys.argv[3])
nt = int(sys.argv[4])
dt = float(sys.argv[5])
ntiskp = int(sys.argv[6])
ts_per_file = int(sys.argv[7])
prefix = sys.argv[8]
output_filename = sys.argv[9]

with open(output_filename, "w") as fp_out:
	fp_out.write('<?xml version="1.0" ?>\n')
	fp_out.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
	fp_out.write('<Xdmf xmlns:xi="http://www.x3.org/2001/XInclude" Version="3.0">\n')
	fp_out.write('\t<Domain>\n')
	fp_out.write('\t\t<Topology name="topo" TopologyType="2DCoRectMesh" Dimensions="%d %d">\n' % (y_dim, x_dim))
	fp_out.write('\t\t</Topology>\n')
	fp_out.write('\t\t<Geometry GeometryType="Origin_DxDy">\n')
	fp_out.write('\t\t\t<DataItem Format="XML" NumberType="Float" Dimensions="2">\n')
	fp_out.write('\t\t\t\t0 0\n')
	fp_out.write('\t\t\t</DataItem>\n')
	fp_out.write('\t\t\t<DataItem Format="XML" NumberType="Float" Dimensions="2">\n')
	fp_out.write('\t\t\t\t%f %f\n' % (grid_dec, grid_dec))
	fp_out.write('\t\t\t</DataItem>\n')
	fp_out.write('\t\t</Geometry>\n')
	fp_out.write('\t\t<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n')
	fp_out.write('\t\t\t<Time TimeType="HyperSlab">\n')
	fp_out.write('\t\t\t\t<DataItem Format="XML" NumberType="Float" Dimensions="3">\n')
	fp_out.write('\t\t\t\t\t%f %f %d\n' % (dt*ntiskp, dt*ntiskp, nt/ntiskp))
	fp_out.write('\t\t\t\t</DataItem>\n')
	fp_out.write('\t\t\t</Time>\n')
	#Section for each timestep
	for i in range(0, nt/ntiskp):
		current_ts = (i+1)*dt*ntiskp
		current_ts_filename = min(((i)/ts_per_file + 1)*ts_per_file*ntiskp, nt)
		offset = 4*x_dim*y_dim*(i%ts_per_file)
		input_filename = "%s%07d" % (prefix, current_ts_filename)
		fp_out.write('\t\t\t<Grid Name="mesh" GridType="Uniform">\n')
		fp_out.write('\t\t\t\t<Topology Reference="/Xdmf/Domain/Topology[1]"/>\n')
		fp_out.write('\t\t\t\t<Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n')
		fp_out.write('\t\t\t\t<Time Value="%f" />\n' % (current_ts))
		fp_out.write('\t\t\t\t<Attribute Name="X Velocity" Center="Node">\n')
		fp_out.write('\t\t\t\t\t<DataItem DataType="Float" Precision="4" Dimensions="%d %d" Format="Binary" Seek="%d" Endian="Little">\n' % (y_dim, x_dim, offset))
		fp_out.write('\t\t\t\t\t\t%s\n' % (input_filename))
		fp_out.write('\t\t\t\t\t</DataItem>\n')
		fp_out.write('\t\t\t\t</Attribute>\n')
		fp_out.write('\t\t\t</Grid>\n')
	fp_out.write('\t\t</Grid>\n')
	fp_out.write('\t</Domain>\n')
	fp_out.write('</Xdmf>\n')
	fp_out.flush()
	fp_out.close()

