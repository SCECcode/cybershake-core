#!/usr/bin/env python3

import sys
import os
import struct
import numpy as np

sys.path.append("/work2/00349/scottcal/frontera/CyberShake/software/DirectSynth/src/duration")

from record_scott import Record

DAMPING = 0.05
FMIN = 0.01
N_HAN_WIN = 13

def populate_waveform_arrays(input_filename):
	with open(input_filename, "rb") as fp_in:
		header_str = fp_in.read(56)
		rv = struct.unpack('i', header_str[32:36])[0]
		dt = struct.unpack("f", header_str[36:40])[0]
		nt = struct.unpack("i", header_str[40:44])[0]
		comp_int = struct.unpack("i", header_str[44:48])[0]
		num_comps = 0
		if (comp_int & 1)==1:
			x_comp = True
			num_comps += 1
		if (comp_int & 2)==2:
			y_comp = True
			num_comps += 1
		if (comp_int & 4)==4:
			z_comp = True
			num_comps += 1
		#Only want horizontal components
		if num_comps==3:
			fp_in.seek(4*nt, os.SEEK_CUR)
		x_str = fp_in.read(4*nt)
		x_data = struct.unpack('%df' % nt, x_str)
		y_str = fp_in.read(4*nt)
		y_data = struct.unpack('%df' % nt, y_str)
		times = [dt*i for i in range(0, nt)]
		return (times, x_data, y_data)


#This function constructs a Record object so that c_calc_rec_durations can be called in a loop over all the frequencies and durations without needing additional record objects.
def construct_rec(seis, times):
	rec = Record(np.array(seis), times, f_min=FMIN, f_max=0.5/(times[1]-times[0])-.5, dt=(times[1]-times[0]), sim_time=times[-1], n_han_win=N_HAN_WIN)
	return rec

#This function returns a single duration value for a given frequency, duration metric, and Record.
def c_calc_rec_durations(freq, dur, rec):
	(xi, xf) = [int(x) for x in dur.split("-")]
	dur_val = round(float(rec.duration(freq, xi=xi, xf=xf, damping=DAMPING)), 2)
	return dur_val

#This function constructs the Record and then calculates the duration.
def c_calc_float_durations(freq, seis, dur, times):
	rec = Record(np.array(seis), times, f_min=FMIN, f_max=0.5/(times[1]-times[0])-.5, dt=(times[1]-times[0]), sim_time=times[-1], n_han_win=N_HAN_WIN)
	(xi, xf) = [int(x) for x in dur.split("-")]
	dur_val = float(rec.duration(freq, xi=xi, xf=xf, damping=DAMPING))
	return dur_val

#This code takes a list of frequencies and a list of durations, calls c_calc_float_durations() to do the calculations, then prints the results.
def calculate_standalone_durations(freqs, durations, time, u, v):
	#To avoid filtering, set fmax to just under the nyquist freq, 0.5/dt
	#u_rec = Record(np.array(u), time, f_min=FMIN, f_max=0.5/(time[1]-time[0])-.01, dt=(time[1]-time[0]), sim_time=time[-1], n_han_win=N_HAN_WIN)
	#v_rec = Record(np.array(v), time, f_min=FMIN, f_max=0.5/(time[1]-time[0])-.01, dt=(time[1]-time[0]), sim_time=time[-1], n_han_win=N_HAN_WIN)
	for f in freqs:
		for d in durations:
			x_dur_val = c_calc_float_durations(f, u, d, time)
			y_dur_val = c_calc_float_durations(f, v, d, time)
			#(xi, xf) = [int(x) for x in d.split("-")]
			#x_dur_val = u_rec.duration(f, xi=xi, xf=xf, damping=DAMPING)
			#y_dur_val = v_rec.duration(f, xi=xi, xf=xf, damping=DAMPING)
			print("Velocity duration (%s) for frequency %f, x-component = %f" % (d, f, x_dur_val))
			print("Velocity duration (%s) for frequency %f, y-component = %f" % (d, f, y_dur_val))

#This is the test created by Camillo for verification.
def default_duration_calcs():
    from scipy import integrate
	# LOADING WAVEFORMS
    path = 'J056_37.77834_-122.44289.txt'
    rec = np.genfromtxt(path)
    time, u, v, w = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3]
    # WAVEFORMS ARE IN ACCELERATION, SO I CONVERT THEM TO VELOCITY TO INITIATE THE MODULE
    vel = integrate.cumtrapz(u, time, initial=0)
    dt = (time[1:] - time[:-1]).mean()

    # SET UP PARAMETERS FOR SDOF SYSTEM: FREQUENCY AND DAMPING
    osc_freqs, damping = 0.1, 0.05
    #osc_freqs = [0.1, 0.2, 0.5, 1.0, 2.0]
    #damping = 0.05
    # LOAD MODULE: THIS MODULE APPLIES A WAVEFORM FILTERING, CONSIDERING THAT IN SYNTHETIC WAVEFORMS, WE HAVE A MAXIMUM FREQUENCY.
    # DT IN THE MODULE IS THE RESAMPLING DT, TO HOMOGENIZE THE TIME STEPS
    fmin, fmax, dt_homogenize = 0.01, 49.5, 0.01
    test_rec = Record(vel, time, f_min=fmin, f_max=fmax, dt=dt_homogenize, sim_time=time[-1], n_han_win=13)
    # SDOF RESPONSE: VECTORS IN CASE YOU WANT TO PLOT THE SDOF RESPONSE
    for o in [0.1, 0.2, 0.5, 1.0, 2.0]:
        osc_freqs = o
        sdof_resp = test_rec.sdof_response(osc_freqs, damping=damping)
        time_sdof = test_rec.time_vector

        # HERE, WE COMPUTE THE DURATION BETWEEN TWO TEST ENERGY LEVELS, Xi AND Xf. FOR D5-95, USE Xi=5 AND Xf=95
        print('Dur: {:.3}'.format(test_rec.duration(osc_freqs, xf=95, xi=5, damping=0.05)))

