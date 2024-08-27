#!/usr/bin/env python3

import sys
import os
import pymysql

if len(sys.argv)<4:
    print("Usage: %s <seeds input file> <erf id> <rupture variation scenario id>" % ss.argv[0])
    sys.exit(1)

seeds_file = sys.argv[1]
erf_id = int(sys.argv[2])
rup_var_scenario_id = int(sys.argv[3])

conn = pymysql.connect(user='cybershk', passwd='***REMOVED***', db='CyberShake', host='moment.usc.edu')
cur = conn.cursor()

with open(seeds_file, "r") as fp_in:
    data = fp_in.readlines()
    fp_in.close()
    for i, line in enumerate(data):
        pieces = line.split()
        src_id = int(pieces[0])
        rup_id = int(pieces[1])
        rv_id = int(pieces[2])
        seed = int(pieces[3])
        query = "insert into Rup_Var_Seeds (ERF_ID, Rup_Var_Scenario_ID, Source_ID, Rupture_ID, Rup_Var_ID, Rup_Var_Seed) values (%d, %d, %d, %d, %d, %d)" % (erf_id, rup_var_scenario_id, src_id, rup_id, rv_id, seed)
        cur.execute(query)
conn.commit()
conn.close()
