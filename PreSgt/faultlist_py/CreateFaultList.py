#!/usr/bin/env python

import os
import sys

import MySQLdb

#server = "focal.usc.edu"
server = "moment.usc.edu"
db = "CyberShake"
user = "cybershk_ro"
passwd = "CyberShake2007"

if len(sys.argv)<5:
	print "Usage: %s <site> <erf_id> <path to rup vars> <output_file> [-rsqsim]" % sys.argv[0]
	sys.exit(1)

site = sys.argv[1]
erf_id = int(sys.argv[2])
rv_path = sys.argv[3]
output_file = sys.argv[4]
#Using RSQSim-formatted rupture geometries (average area, only num points, so only 4 rows in header)
header_rows = 6

if len(sys.argv)==6:
	if sys.argv[5]=="-rsqsim":
		header_rows = 4

conn = MySQLdb.connect(host=server, user=user, passwd=passwd, db=db)
cur = conn.cursor()
query = 'select R.Source_ID, R.Rupture_ID from CyberShake_Site_Ruptures R, CyberShake_Sites S where S.CS_Short_Name="%s" and R.CS_Site_ID=S.CS_Site_ID and R.ERF_ID=%d order by R.Source_ID, R.Rupture_ID' % (site, erf_id)
cur.execute(query)
res = cur.fetchall()
if len(res)==0:
	print "No ruptures found for site %s." % site
	sys.exit(2)

fp_out = open(output_file, "w")
count = 0
for r in res:
	sourceID = int(r[0])
	ruptureID = int(r[1])
	fp_out.write("%s/%d/%d/%d_%d.txt nheader=%d latfirst=1\n" % (rv_path, sourceID, ruptureID, sourceID, ruptureID, header_rows))
	count += 1
	if (count%100==0):
		print "Processed %d ruptures." % count
fp_out.flush()
fp_out.close()

cur.close()
conn.close()
