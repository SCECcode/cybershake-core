#!/usr/bin/env python3

'''This script takes a broadband CyberShake run ID, an LF run ID, and a BB velocity info file.  It then inserts Vs30 from the velocity info file into the database for the broadband run as Model_Vs30, and copies Z1.0 and Z2.5 over from the LF run, if they exist.'''

import sys
import os
import pymysql
import argparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

parser = argparse.ArgumentParser(description = "This script collects Vs30 from a Broadband CyberShake velocity info file, and inserts it into the database for the broadband run as 'Model_Vs30'.  It also copies Z1.0 and Z2.5 over from the LF run to the BB run.")
parser.add_argument("-vi", "--velocity_info", help="Path to BB CS velocity info file (required)", type=str)
parser.add_argument("-bbid", "--broadband_run_id", help="Run ID of the broadband CyberShake run (required)", type=int)
parser.add_argument("-lfid", "--lowfreq_run_id", help="Run ID of the low-frequency CyberShake run (required)", type=int)

args = parser.parse_args()
if args.velocity_info is None or args.broadband_run_id is None or args.lowfreq_run_id is None:
    print("Missing required arguments.")
    parser.print_help()
    sys.exit(1)

vel_info_file = args.velocity_info
bbid = args.broadband_run_id
lfid = args.lowfreq_run_id

#Get Vs30 from vel info file
with open(vel_info_file, 'r') as fp_in:
    data = fp_in.readlines()
    vs30 = -1.0
    for line in data:
        pieces = line.split("=")
        if pieces[0].strip()=='Vs30':
            vs30 = float(pieces[1])
            break
    fp_in.close()
if vs30<0:
    print("Vs30 not found in velocity info file %s, aborting." % vel_info_file)
    sys.exit(2)

#Connect to DB
db_file = config.getProperty("DB_WR_FILE")
with open(db_file, "r") as fp_in:
    username = fp_in.readline().strip()
    password = fp_in.readline().strip()
    fp_in.close()

conn = pymysql.connect(host='moment.usc.edu', db='CyberShake', passwd=password, user=username)
cur = conn.cursor()

#Get Z1.0, Z2.5 from LF run
z1_0 = None
z2_5 = None
query = 'select Z1_0, Z2_5 from CyberShake_Runs where Run_ID=%d' % (lfid)
cur.execute(query)
res = cur.fetchone()
if res[0] is not None:
    z1_0 = float(res[0])
if res[1] is not None:
    z2_5 = float(res[1])

if z1_0 is None:
    z1_0_str = 'NULL'
else:
    z1_0_str = str(z1_0)
if z2_5 is None:
    z2_5_str = 'NULL'
else:
    z2_5_str = str(z2_5)

query = 'update CyberShake_Runs set Model_Vs30=%f, Z1_0=%s, Z2_5=%s where Run_ID=%d' % (vs30, z1_0_str, z2_5_str, bbid)
print(query)
cur.execute(query)
conn.commit()
conn.close()
