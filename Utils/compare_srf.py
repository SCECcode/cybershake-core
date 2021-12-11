#!/usr/bin/env python

import sys
import os
import math

def cmp_srf(filename1, filename2, tolerance=0.0011):
	'''
	Compare two SRF files.  A tolerance is accepted for floating-point
	values.  If any point parameters differ by more than epsilon
	but less than the tolerance, let it go, but don't compare the slip
	values for that point (since they may differ by quite a bit).
	Count up the number of times it happens;  if it's more than 1% of the
	number of points reject it.
	'''

	fp1 = open(filename1, 'r')
	fp2 = open(filename2, 'r')
	data1 = fp1.readlines()
	data2 = fp2.readlines()
	fp1.close()
	fp2.close()

	returncode = 0

	if not data1[0]==data2[0]:
		#incompatible versions
		print "SRF versions don't match."
		return 1

	i1 = 1
	i2 = 1
	line1 = data1[i1]
	line2 = data2[i2]
	pieces1 = line1.split()
	#first PLANE section
	if pieces1[0]=="PLANE":
		if not line1==line2:		
			print "Line 0: mismatched planes."
			return 1
		num_planes = int(pieces1[1])
		#for each plane
		for j in range(0, num_planes):
			#for each line of plane data
			for k in range(0, 2):
				if not data1[i1+2*j+k]==data2[i2+2*j+k]:
					print "Line %d: plane data doesn't match." % (i1+2*j+k)
					return 1
		i1 = i1+1+2*num_planes
		i2 = i2+1+2*num_planes
	else:
		print "Malformed SRF file;  no PLANEs"
		return -1
	#then POINTS
	NT_INDEX = 2
	NT_TOLERANCE = 1
	LL_TOLERANCE = 0.0001
	PP_TOLERANCE = 0.031
	RAKE_TOLERANCE = 1
	SKIP_TOLERANCE = 0.0001
	max_diff = 0.0
	max_percent_diff = 0.0
	line1 = data1[i1]
	line2 = data2[i2]
	pieces1 = line1.split()
	if pieces1[0]=="POINTS":
		if not line1==line2:
			print "Line %d: points line doesn't match." % (i1)
			return 1
		num_points = int(pieces1[1])
		points_skipped = 0
		for j in range(0, num_points):
			#if any of the params are different, skip the slip values
			skip_slips = False
			#each point has 2 lines + ceil(NT1/6) lines
			#will permit tolerance in lat/lon points
			#and tolerance in NT values
			line1 = data1[i1+1]
			line2 = data2[i2+1]
			pieces1 = line1.split()
			pieces2 = line2.split()
			#check line length
			if not len(pieces1)==len(pieces2):
				print "Line %d (point %d):  files have different number of point parameters." % (i1+1, j)
				print line1, line2
				return 1
			#check lons
			lon1 = float(pieces1[0])
			lon2 = float(pieces2[0])
			if math.fabs(lon1-lon2)>LL_TOLERANCE:
				print "Line %d/%d (point %d):  longitudes %f and %f differ by more than the accepted tolerance %f." % (i1+1, i2+1, j, lon1, lon2, LL_TOLERANCE)
				return 2
			#check lats
                        lat1 = float(pieces1[1])
                        lat2 = float(pieces2[1])
                        if math.fabs(lat1-lat2)>LL_TOLERANCE:
                                print "Line %d/%d (point %d):  latitudes %f and %f differ by more than the accepted tolerance %f." % (i1+1, i2+1, j, lat1, lat2, LL_TOLERANCE)
                                return 3
			#check rest of line
			for k in range(2, len(pieces1)):
				if math.fabs(float(pieces1[k])-float(pieces2[k]))>PP_TOLERANCE:
					print "Line %d/%d (point %d):  point parameters in field %d disagree." % (i1+1, i2+1, j, k+1)
					print line1, line2
					return 4
				if math.fabs(float(pieces1[k])-float(pieces2[k]))>SKIP_TOLERANCE:
					#don't compare the slip values, they'll be different
					skip_slips = True
			#second line
			line1 = data1[i1+2]
			line2 = data2[i2+2]
                        pieces1 = line1.split()
                        pieces2 = line2.split()
			#check line length
                        if not len(pieces1)==len(pieces2):
                                print "Line %d/%d (point %d):  files have different number of point parameters." % (i1+2, i2+2, j)
                                print line1, line2
                                return 1
			#check rake
			r1 = int(pieces1[0])
			r2 = int(pieces2[0])
			if r1!=r2:
				skip_slips = True
			if abs(r1-r2)>RAKE_TOLERANCE:
				print "Line %d/%d (point %d):  rakes %d and %d differ by more than the accepted tolerance of 1." % (i1+2, i2+2, j, r1, r2)
                                print line1, line2
                                return 4
			for k in range(1, NT_INDEX):
				if math.fabs(float(pieces1[k])-float(pieces2[k]))>PP_TOLERANCE:
                                        print "Line %d/%d (point %d):  point parameters in field %d disagree." % (i1+2, i2+2, j, k+1)
                                        print line1, line2
                                        return 4
				if math.fabs(float(pieces1[k])-float(pieces2[k]))>SKIP_TOLERANCE:
					#ok, but skip comparisons
					skip_slips = True
			#compare NTs
			nt1 = int(pieces1[NT_INDEX])
			nt2 = int(pieces2[NT_INDEX])
			if not nt1==nt2:
				skip_slips = True
				if abs(nt1-nt2)>NT_TOLERANCE:
					print "Line %d/%d (point %d):  NT values %d and %d differ by %f, more than the accepted tolerance %d." % (i1+2, i2+2, j, nt1, nt2, abs(nt1-nt2), NT_TOLERANCE)
					return 5
			if not int(pieces1[4])==0 or not int(pieces1[6])==0:
				print "Line %d: SRF has NT2 or NT3, need to alter parser." % (i1+2)
				sys.exit(0)
			if skip_slips==False:
				num_rows = int(math.ceil(nt1/6.0))
				for k in range(0, num_rows):
					line1 = data1[i1+3+k]
					line2 = data2[i2+3+k]
					pieces1 = line1.split()
					pieces2 = line2.split()
					if not len(pieces1)==len(pieces2):
						print "Line %d/%d: mismatch in entries in line." % (i1+3+k, i2+3+k)
						continue
					for p in range(0, len(pieces1)):
						p1 = float(pieces1[p])
						p2 = float(pieces2[p])
						if p1<1.0 or p2<1.0:
							if math.fabs(p1-p2)>tolerance:
								print "Line %d/%d (point %d): %f and %f differ by more than the accepted tolerance %f (%f)." % (i1+3+k, i2+3+k, j, float(pieces1[p]), float(pieces2[p]), tolerance, math.fabs(p1-p2))
								returncode = 1
							if math.fabs(p1-p2)>max_diff:
								max_diff = math.fabs(p1-p2)
						else:
							if math.fabs(p1-p2)/p1>tolerance:
								print "Line %d/%d (point %d): %f and %f differ by more than the accepted tolerance %f%% (%f%%)." % (i1+3+k, i2+3+k, j, float(pieces1[p]), float(pieces2[p]), tolerance*100.0, math.fabs(p1-p2)/p1*100.0)
								returncode = 1
							if math.fabs(p1-p2)/p1>max_percent_diff:
								max_percent_diff = math.fabs(p1-p2)/p1
				i1 += 2+num_rows
				i2 += 2+num_rows
			else:
				points_skipped += 1
				i1 += 2+int(math.ceil(nt1/6.0))
				i2 += 2+int(math.ceil(nt2/6.0))
	else:
		print "Malformed SRF file;  no POINTS."
		return -2
	if points_skipped > num_points/50:
		print "Too many points with different parameters;  of %d total points %d had different parameters." % (num_points, points_skipped)
		returncode = 2

	print "Max diff: %f, max percent diff: %f%%" % (max_diff, max_percent_diff*100.0)
	return returncode
		

def cmp_resid(filename1, filename2, tolerance=0.0015):
	'''
	Format is long list of heaters, then for each station
	<eq> <mag> <stat name> <lon> <lat> <stat_seq_no> <vs30> <close_dist> <Xcos> <Ycos> <T_min> <T_max> <comp> <period1> ...
	'''
	fp1 = open(filename1, 'r')
	fp2 = open(filename2, 'r')
	data1 = fp1.readlines()
	data2 = fp2.readlines()
	fp1.close()
	fp2.close()
	
	returncode = 0

	i = 0
	
	line1 = data1[i]
	line2 = data2[i]
	if line1!=line2:
		print "Discrepancy in header lines."
		returncode = -1
		return returncode
	for i in range(1, len(data1)):
		pieces1 = data1[i].split()
		pieces2 = data2[i].split()
		for j in range(0, 13):
			if pieces1[j]!=pieces2[j]:
				print "Line %d: %s and %s don't agree." % ((i+1), pieces1[j], pieces2[j])
		for j in range(13, len(pieces1)):
			f1 = float(pieces1[j])
			f2 = float(pieces2[j])
			if math.fabs(f1)<1.0:
				if math.fabs(f1-f2)>tolerance:
					print "Line %d: %f and %f differ by more than %f tolerance." % ((i+1), f1, f2, tolerance)
					returncode = 1
			else:
				if math.fabs(f1-f2)/f1>tolerance:
					print "Line %d:  %f and %f differ by more than %f%% tolerance." % ((i+1), f1, f2, tolerance*100.0)
                                        returncode = 1
	return returncode

	

def cmp_bbp(filename1, filename2, tolerance=0.0015):
	fp1 = open(filename1, 'r')
	fp2 = open(filename2, 'r')
	data1 = fp1.readlines()
	data2 = fp2.readlines()
	fp1.close()
	fp2.close()

	returncode = 0

	#need to skip comments in both files
	file1_offset = 0
	file2_offset = 0
	for line1 in data1:
		if line1.strip().startswith("#") or line1.strip().startswith("%"):
			file1_offset += 1
		else:
			break
        for line2 in data2:
                if line2.strip().startswith("#") or line2.strip().startswith("%"):
                        file2_offset += 1
                else:
                        break


	for i in range(0, len(data1)-file1_offset):
		line1 = data1[i+file1_offset]
		pieces1 = line1.split()
		line2 = data2[i+file2_offset]
		pieces2 = line2.split()
		if not float(pieces1[0])==float(pieces2[0]): #timestamps really should be equal
			print "Timestamps %f in file %s and %f in file %s aren't equal." % (float(pieces1[0]), filename1, float(pieces2[0]), filename2)
			returncode = -1
			return returncode
		for j in range(1, 4):
			f1 = float(pieces1[j])
			f2 = float(pieces2[j])
			if math.fabs(f1)<1.0 or math.fabs(f2)<1.0:
                        	if math.fabs(f1-f2)>tolerance:
                                	print "Line %d:  %f and %f differ by more than %f tolerance." % ((i+1), f1, f2, tolerance)
                                	returncode = 1
			else:
				if math.fabs(f1-f2)/f1>tolerance:
					print "Line %d:  %f and %f differ by more than %f%% tolerance." % ((i+1), f1, f2, tolerance*100.0)
					returncode = 1
		if i>1000:
			return returncode
	return returncode


def cmp_bias(filename1, filename2, tolerance=0.0015):
        fp1 = open(filename1, 'r')
        fp2 = open(filename2, 'r')
        data1 = fp1.readlines()
        data2 = fp2.readlines()
        fp1.close()
        fp2.close()

        returncode = 0

	for i in range(0, len(data1)):
		pieces1 = data1[i].split()
		pieces2 = data2[i].split()
		f1 = float(pieces1[1])
		f2 = float(pieces2[1])

		if pieces1[0]!=pieces2[0]:
			print "Line %d:  periods %f and %f don't agree." % ((i+1), float(pieces1[0]), float(pieces2[0]))
			returncode = 1
		if math.fabs(f1)<1.0:
			if math.fabs(f1-f2)>tolerance:
	                        print "Line %d:  %f and %f differ by more than %f tolerance." % ((i+1), f1, f2, tolerance)
                                returncode = 1
                else:
                        if math.fabs(f1-f2)/f1>tolerance:
                                print "Line %d:  %f and %f differ by more than %f%% tolerance." % ((i+1), f1, f2, tolerance*100.0)
                                returncode = 1
	return returncode


cmp_srf(sys.argv[1], sys.argv[2], tolerance=0.00015)
