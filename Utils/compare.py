#!/usr/bin/env python

import os
import sys
import struct

if len(sys.argv)<3:
	print "Usage: %s <file1> <file2>" % sys.argv[0]
	sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]

fp1 = open(file1, "rb")
fp2 = open(file2, "rb")

data1 = fp1.read(4)
data2 = fp2.read(4)

totDiff = 0.0
totPercentDiff = 0.0
biggestDiff = 0.0
biggestDiffPercent = 0.0
biggestDiffIndex = -1
biggestPercentDiff = 0.0
biggestPercentDiffDiff = 0.0
biggestPercentDiffIndex = -1

diffLargerThanOne = []
indexDiffLargerThanOne = []

counter = 0

while not data1=="":
	float1 = float(struct.unpack("f", data1)[0])
        float2 = float(struct.unpack("f", data2)[0])
	
	diff = abs(float2-float1)
        if diff>biggestDiff:
                biggestDiff = diff
                biggestDiffPercent = percentDiff
                biggestDiffIndex = counter
	if (float1!=0.0):
		percentDiff = diff/float1 * 100.0
		if percentDiff>biggestPercentDiff:
			biggestPercentDiff = percentDiff
			biggestPercentDiffDiff = diff
			biggestPercentDiffIndex = counter
	if diff>10.0:
		print "Diff of %f in index %d." % (diff, counter)
		sys.exit(1)
		diffLargerThanOne.append(diff)
		indexDiffLargerThanOne.append(counter)
	totDiff += diff
	totPercentDiff += percentDiff
	counter += 1
	if counter % 1000000 == 0:
		print "%d entries." % counter
	data1 = fp1.read(4)
	data2 = fp2.read(4)

fp1.close()
fp2.close()

print "Average diff: %f." % (totDiff/counter)
print "Average percent diff: %f." % (totPercentDiff/counter)
print "Max diff: %f, with percentage %f, index %d." % (biggestDiff, biggestDiffPercent, biggestDiffIndex)
print "Max percent diff: %f, with diff %f, index %d." % (biggestPercentDiff, biggestPercentDiffDiff, biggestPercentDiffIndex)

#fp_out = open("ucvm_diffs.txt", "w")
#for i in range(0, len(diffLargerThanOne)):
#	fp_out.write("index %d, diff %f\n" % (indexDiffLargerThanOne[i], diffLargerThanOne[i]))
#fp_out.flush()
#fp_out.close()

