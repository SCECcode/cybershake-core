#!/usr/bin/env python

import os
import sys
import optparse
import pymysql
import sqlite3

#server = "focal.usc.edu"
server = "moment.usc.edu"
db = "CyberShake"
user = "cybershk_ro"
passwd = "CyberShake2007"

parser = optparse.OptionParser()
parser.add_option("--site", dest="site", action="store", help="Site name")
parser.add_option("--erf_id", dest="erf_id", action="store", type="int", help="ERF ID")
parser.add_option("--rup_path", dest="rup_path", action="store", help="Path to root directory containing rupture files.")
parser.add_option("--output", dest="output_file", action="store", help="Path to output file.")
parser.add_option("--server", dest="server", action="store", help="Path to server (default is %s).  Can specify SQLite file with sqlite://<file>")
parser.add_option("--rsqsim", dest="rsqsim", action="store_true", default=False, help="Assumes RSQSim-formatted rupture geometry files.  Default is false.")

(option, args) = parser.parse_args()

site = option.site
erf_id = option.erf_id
rv_path = option.rup_path
output_file = option.output_file
#Using RSQSim-formatted rupture geometries (average area, only num points, so only 4 rows in header)
header_rows = 6

if option.rsqsim==True:
	header_rows = 4

if (site==None or erf_id==None or rv_path==None or output_file==None):
	print("site, erf_id, rup_path, and output are all required.")
	parser.print_help()
	sys.exit(1)

if option.server is not None:
	server = option.server

if server[0:9]=="sqlite://":
	conn = sqlite3.connect(server[9:])
else:
	conn = pymysql.connect(host=server, user=user, passwd=passwd, db=db)
cur = conn.cursor()
query = 'select R.Source_ID, R.Rupture_ID from CyberShake_Site_Ruptures R, CyberShake_Sites S where S.CS_Short_Name="%s" and R.CS_Site_ID=S.CS_Site_ID and R.ERF_ID=%d order by R.Source_ID, R.Rupture_ID' % (site, erf_id)
cur.execute(query)
res = cur.fetchall()
if len(res)==0:
	print("No ruptures found for site %s." % site)
	sys.exit(2)

fp_out = open(output_file, "w")
count = 0
for r in res:
	sourceID = int(r[0])
	ruptureID = int(r[1])
	fp_out.write("%s/%d/%d/%d_%d.txt nheader=%d latfirst=1\n" % (rv_path, sourceID, ruptureID, sourceID, ruptureID, header_rows))
	count += 1
	if (count%100==0):
		print("Processed %d ruptures." % count)
fp_out.flush()
fp_out.close()

cur.close()
conn.close()
