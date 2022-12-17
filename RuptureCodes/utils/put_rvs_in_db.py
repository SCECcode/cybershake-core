#!/usr/bin/env python3

import sys
import os
import pymysql

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))
sys.path.append(path_add)

import config

if len(sys.argv)<4:
    print("Usage: %s <erf id> <rup var scenario id> <rv list file>" % sys.argv[0])
    sys.exit(1)

erf_id = int(sys.argv[1])
rup_var_scenario_id = int(sys.argv[2])
rv_list_file = sys.argv[3]

db_pass_file = config.getProperty("DB_WR_FILE")
with open(db_pass_file, "r") as fp_in:
    user = fp_in.readline().strip() 
    passwd = fp_in.readline().strip()
    fp_in.close()

conn = pymysql.connect(user=user, passwd=passwd, host='moment.usc.edu', db='CyberShake')
cur = conn.cursor()

with open(rv_list_file, "r") as fp_in:
    data = fp_in.readlines()
    fp_in.close()
    for i in range(0, len(data), 2):
        if (i%200==0):
            print("Processing rupture %d of %d." % (i//2, len(data)//2))
        (src_id, rup_id) = data[i].split()
        num_rvs = int(data[i+1])
        for j in range(0, num_rvs):
            query = 'insert into Rupture_Variations(Rup_Var_ID, Rup_Var_Scenario_ID, ERF_ID, Source_ID, Rupture_ID, Rup_Var_LFN) values (%d, %d, %d, %s, %s, "e%d_rv%d_%s_%s.txt.variation-r000%03d")' % (j, rup_var_scenario_id, erf_id, src_id, rup_id, erf_id, rup_var_scenario_id, src_id, rup_id, j)
            cur.execute(query)
        conn.commit()
cur.close()
conn.close()
