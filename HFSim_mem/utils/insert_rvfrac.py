#!/usr/bin/env python3

import sys
import os
import pymysql

if len(sys.argv)<4:
    print("Usage: %s <rvfrac values file> <erf id> <rup var scenario id>" % sys.argv[0])
    sys.exit(1)

conn = pymysql.connect(user='cybershk', passwd='***REMOVED***', db='CyberShake', host='moment.usc.edu')
cur = conn.cursor()

rvfrac_file = sys.argv[1]
erf_id = int(sys.argv[2])
rup_var_scenario_id = int(sys.argv[3])

#Lines look like 0 0 0 0.843884
with open(rvfrac_file, "r") as fp_in:
    data = fp_in.readlines()
    fp_in.close()
    for i, line in enumerate(data):
        #if (i%100==0):
            #print("Line %d of %d" % (i, len(data)))
        pieces = line.split()
        src_id = int(pieces[0])
        rup_id = int(pieces[1])
        rv_id = int(pieces[2])
        rvfrac = float(pieces[3])
        query = "update Rupture_Variations set rvfrac=%f where Rup_Var_ID=%d and Rup_Var_Scenario_ID=%d and ERF_ID=%d and Source_ID=%d and Rupture_ID=%d" % (rvfrac, rv_id, rup_var_scenario_id, erf_id, src_id, rup_id)
        cur.execute(query)
conn.commit()
conn.close()
