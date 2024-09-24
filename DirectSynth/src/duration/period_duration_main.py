#!/usr/bin/env python3

import sys
import numpy as np
from scipy import integrate
import struct

sys.path.append("/work2/00349/scottcal/frontera/CyberShake/software/DirectSynth/src/duration")

from record_scott import Record
import period_duration

def run_main(argv=sys.argv):
	#period_duration.default_duration_calcs()
	#sys.exit(0)
	if len(sys.argv)<4:
		print("Usage: %s <single-RV .grm file> <comma-separated periods to calculate duration> <comma-separated durations>" % (sys.argv[0]))
		print("Example: %s Seismogram_TEST_0_0.grm 2,3,5,7.5,10 5-75,5-95,20-80" % sys.argv[0])
		sys.exit(1)
	input_filename = sys.argv[1]
	periods = [float(p) for p in sys.argv[2].split(",")]
	osc_freqs = [1.0/p for p in periods]
	durations = sys.argv[3].split(",")
	(time, u, v) = period_duration.populate_waveform_arrays(input_filename)
	period_duration.calculate_standalone_durations(osc_freqs, durations, time, u, v)

	u_rec = period_duration.construct_rec(u, time)
	v_rec = period_duration.construct_rec(v, time)
	for f in osc_freqs:
		print(period_duration.c_calc_rec_durations(f, durations[0], u_rec))
		print(period_duration.c_calc_rec_durations(f, durations[0], v_rec))

if __name__=='__main__':
	run_main()
