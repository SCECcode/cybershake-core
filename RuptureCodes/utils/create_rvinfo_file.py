#!/usr/bin/env python3

import sys
import os
import pymysql

if len(sys.argv)<5:
    print("Usage: %s <site> <erf ID> <rup var scenario id> <outfile>" % (sys.argv[0]))
    sys.exit(1)

site = sys.argv[1]
erf_id = int(sys.argv[2])
rup_var_scenario_id = int(sys.argv[3])
outfile = sys.argv[4]

conn = pymysql.connect(host='moment.usc.edu', user='cybershk_ro', passwd='CyberShake2007', db='CyberShake')
cur = conn.cursor()

query = "select V.Source_ID, V.Rupture_ID, V.Rup_Var_ID, V.rvfrac, D.Rup_Var_Seed " \
"from Rupture_Variations V, Rup_Var_Seeds D, CyberShake_Site_Ruptures SR, CyberShake_Sites S " \
"where S.CS_Short_Name='%s' " \
"and S.CS_Site_ID=SR.CS_Site_ID " \
"and SR.ERF_ID=%d " \
"and SR.ERF_ID=V.ERF_ID " \
"and SR.Source_ID=V.Source_ID " \
"and SR.Rupture_ID=V.Rupture_ID " \
"and SR.Cutoff_Dist=200.0 " \
"and V.Rup_Var_Scenario_ID=%d " \
"and D.Rup_Var_Scenario_ID=V.Rup_Var_Scenario_ID " \
"and D.ERF_ID=SR.ERF_ID " \
"and D.Source_ID=V.Source_ID " \
"and D.Rupture_ID=V.Rupture_ID " \
"and D.Rup_Var_ID=V.Rup_Var_ID " \
"order by D.Source_ID asc, D.Rupture_ID asc, D.Rup_Var_ID asc" % (site, erf_id, rup_var_scenario_id)
print(query)

cur.execute(query)
res = cur.fetchall()

with open(outfile, 'w') as fp_out:
    fp_out.write("%d\n" % len(res))
    for r in res:
        fp_out.write("%d %d %d %f %d\n" % (int(r[0]), int(r[1]), int(r[2]), float(r[3]), int(r[4])))
    fp_out.flush()
    fp_out.close()

conn.close()

