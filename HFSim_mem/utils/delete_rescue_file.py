#!/usr/bin/env python3

import sys
import os

if len(sys.argv)<4:
    print("Usage: %s <merge pmc rescue file> <CSV of src,rup to delete> <output rescue file>" % sys.argv[0])
    sys.exit(1)

rescue_in = sys.argv[1]
rupture_file = sys.argv[2]
rescue_out = sys.argv[3]

tasks = []

with open(rescue_in, 'r') as fp_in:
    data = fp_in.readlines()
    for line in data:
        tasks.append(line.strip())
    fp_in.close()

rup_strings = []

with open(rupture_file, 'r') as fp_in:
    data = fp_in.readlines()
    for line in data:
        (src_id, rup_id) = line.strip().split(",")
        rup_str = "_%s_%s_" % (src_id, rup_id)
        rup_strings.append(rup_str)
    fp_in.close()

with open(rescue_out, 'w') as fp_out:
    for task in tasks:
        found = False
        for rup_str in rup_strings:
            if task.find(rup_str)>=0:
                found = True
                break
        if found==True:
            continue
        else:
            fp_out.write("%s\n" % task)
    fp_out.flush()
    fp_out.close()

