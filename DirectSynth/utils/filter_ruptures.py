#!/usr/bin/env python3

import sys
import os

if len(sys.argv)<6:
    print("Usage: %s <input rupture list> <input rvfrac list> <CSV with ruptures to include> <output rupture list> <output rvfrac list>" % sys.argv[0])
    sys.exit(1)

input_rup_file = sys.argv[1]
input_rvfrac_file = sys.argv[2]
rup_to_include_file = sys.argv[3]
output_rup_file = sys.argv[4]
output_rvfrac_file= sys.argv[5]

rups_to_include = []

with open(rup_to_include_file, "r") as fp_in:
    data = fp_in.readlines()
    for line in data:
        rups_to_include.append(line.split(',')[0:2])
    fp_in.close()

with open(input_rup_file, "r") as fp_in:
    #Lines look like
    #e36_rv10_128_1296.txt 822 1 120 2740 8.15
    data = fp_in.readlines()
    output_data = []
    for i,line in enumerate(data[1:]):
        if (i%1000==0):
            print("%d of %d" % (i, len(data)))
        pieces = line.split()[0].split("_")
        for (s, r) in rups_to_include:
            if pieces[2]==s and pieces[3].split(".")[0]==r:
                output_data.append(line)
                break
    fp_in.close()

with open(output_rup_file, 'w') as fp_out:
    fp_out.write("%d\n" % len(output_data))
    for line in output_data:
        fp_out.write(line)
    fp_out.flush()
    fp_out.close()

sys.exit(0)
output_data.clear()

with open(input_rvfrac_file, 'r') as fp_in:
    #Lines look like
    #6 0 0 0.804553 60000
    data = fp_in.readlines()
    for i,line in enumerate(data[1:]):
        if (i%1000==0):
            print("%d of %d" % (i, len(data)))
        pieces = line.split()
        for (s, r) in rups_to_include:
            if pieces[0]==s and pieces[1]==r:
                output_data.append(line)
                break
    fp_in.close()

with open(output_rvfrac_file, 'w') as fp_out:
    fp_out.write("%d\n" % len(output_data))
    for line in output_data:
        fp_out.write(line)
    fp_out.flush()
    fp_out.close()


