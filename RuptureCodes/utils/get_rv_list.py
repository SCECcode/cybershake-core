#!/usr/bin/env python3

import sys
import os
import re

if len(sys.argv)<3:
    print("Usage: %s <rupture root dir> <output file>" % sys.argv[0])
    sys.exit(1)

get_num_rvs_path = "/gpfs/alpine/proj-shared/geo112/CyberShake/software/RuptureCodes/RupGen-api-5.5.2/utils/get_num_rvs"
ERF_ID=36

root_dir = sys.argv[1]
output_file = sys.argv[2]

regex = "\d+\d+\.txt$"
p = re.compile(regex)

finished = dict()

with open(output_file, "r") as fp_in:
    data = fp_in.readlines()
    fp_in.close()
    for i in range(0, len(data), 2):
        finished[data[i].strip()] = 1

for source_dir in os.listdir(root_dir):
    src_path = os.path.join(root_dir, source_dir)
    if os.path.isdir(src_path):
        print("Processing src dir %s" % source_dir)
        src_id = int(source_dir)
        for rup_dir in os.listdir(src_path):
            rup_path = os.path.join(src_path, rup_dir)
            if os.path.isdir(rup_path):
                rup_id = int(rup_dir)
                if "%d %d" % (src_id, rup_id) in finished:
                    print("Src %d, rupture %d already done" % (src_id, rup_id))
                    continue
                cmd = "echo %d %d >> %s; %s 36 %d %d | tail -n 1 >> %s" % (src_id, rup_id, output_file, get_num_rvs_path, src_id, rup_id, output_file)
                os.system(cmd)

