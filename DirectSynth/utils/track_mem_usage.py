#!/usr/bin/env python3

import sys
import os

if len(sys.argv)<2:
	print("Usage: %s <direct synth log file>" % sys.argv[0])
	sys.exit(1)

logfile = sys.argv[1]
mem_dict = dict()
num_allocs = 0
num_reallocs = 0
num_frees = 0
with open(logfile, 'r') as fp_in:
	for line in fp_in.readlines():
		if line.find('Allocated')>-1:
			pieces = line.split()
			size = int(pieces[4])
			ptr = pieces[8][:-1]
			mem_dict[ptr] = size
			#print("Allocated %d to ptr %s." % (size, ptr))
		elif line.find('Reallocating')>-1:
			pieces = line.split()
			size = int(pieces[4])
			old_ptr = pieces[9]
			new_ptr = pieces[11][:-1]
			if new_ptr=="175890928":
				print("Reallocated %d from ptr %s to ptr %s." % (size, old_ptr, new_ptr))
			if old_ptr is not '0':
				mem_dict.pop(old_ptr)
			mem_dict[new_ptr] = size
		elif line.find('Freeing')>-1:
			pieces = line.split()
			ptr = pieces[5][:-1]
			#if ptr=="175890928":
			#	print("Freeing ptr %s." % ptr)
			if ptr not in mem_dict:
				print("Can't find pointer %s in mem_dict to free." % ptr)
				continue
				#sys.exit(2)
			mem_dict.pop(ptr)
			#print("Freed ptr %s." % ptr)
	fp_in.close()

#Print everything that's left
print("Remaining allocated memory:")
total = 0
for ptr in mem_dict:
	print("%d bytes at ptr %s" % (mem_dict[ptr], ptr))
	total += mem_dict[ptr]
print("Total: %d MB" % (total/1048576))

