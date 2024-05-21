#!/usr/bin/env python3

'''Create a directory structure for an ERF ID with symlinks which is a subset of another ERF ID.  Both ERF IDs must exist in the database.'''

import sys
import os
import pymysql

if len(sys.argv)<4:
	print("Usage: %s <ERF ID to extract from> <subset ERF ID> <root rupture directory>" % sys.argv[0])
	sys.exit(1)

super_erf_id = int(sys.argv[1])
subset_erf_id = int(sys.argv[2])
root_dir = sys.argv[3]

conn = pymysql.connect(host='moment.usc.edu', user='cybershk_ro', passwd='CyberShake2007', db='CyberShake')
cur = conn.cursor()

query = 'select Source_ID, Rupture_ID from Ruptures where ERF_ID=%d order by Source_ID asc, Rupture_ID asc' % subset_erf_id
cur.execute(query)
res = cur.fetchall()

for r in res:
	(src_id, rup_id) = r
	dst_dir = os.path.join(root_dir, "Ruptures_erf%d" % subset_erf_id, str(src_id), str(rup_id))
	if not os.path.exists(dst_dir):
		os.makedirs(dst_dir)
	sym_src = os.path.join(root_dir, "Ruptures_erf%d" % super_erf_id, str(src_id), str(rup_id), "%d_%d.txt" % (src_id, rup_id))
	sym_dst = os.path.join(dst_dir, "%d_%d.txt" % (src_id, rup_id))
	print("Create symlink %s -> %s" % (sym_dst, sym_src))
	os.symlink(sym_src, sym_dst)

cur.close()
conn.close()
