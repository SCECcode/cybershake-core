#!/usr/bin/env python3

import sys
import os
import pymysql

if len(sys.argv)<4:
    print("Usage: %s <ERF ID> <rupture variation scenario ID> <output file>" % sys.argv[0])
    sys.exit(1)

erf_id = int(sys.argv[1])
rup_var_scenario_id = int(sys.argv[2])
output_file = sys.argv[3]

conn = pymysql.connect(host='moment.usc.edu', user='cybershk_ro', passwd='CyberShake2007', db='CyberShake')
cur = conn.cursor()

query = 'select Source_ID, Rupture_ID, Rup_Var_ID from Rupture_Variations where ERF_ID=%d and Rup_Var_Scenario_ID=%d order by Source_ID asc, Rupture_ID asc, Rup_Var_ID asc' % (erf_id, rup_var_scenario_id)
cur.execute(query)
res = cur.fetchall()

with open(output_file, "w") as fp_out:
    for r in res:
        fp_out.write("%d %d %d\n" % (int(r[0]), int(r[1]), int(r[2])))
    fp_out.flush()
    fp_out.close()
conn.close()
