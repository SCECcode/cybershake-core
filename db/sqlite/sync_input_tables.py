#!/usr/bin/env python3

'''Note: on Frontier this must be run from a head node, because the compute nodes don't have external network access.'''
'''Can process about 300 MB/sec, so plan on 30-45 min to regenerate all tables.'''

import sys
import os
import pymysql
import sqlite3
import argparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(os.path.dirname(full_path)))

sys.path.append(path_add)

import config

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--remote-db', dest="remote_db", action="store", required=True, help="Remote MariaDB or MySQL database")
parser.add_argument('-l', '--local-db', dest="local_db", action="store", required=True, help="Local SQLite DB")
parser.add_argument('-f', '--force', dest='force', action="store_true", required=False, help="Force all local tables to be updated.")

args = parser.parse_args()
remote_db = args.remote_db
local_db = args.local_db
force_update = args.force
if force_update is None:
	force_update = False

my2lite_path = os.path.join(config.getProperty('CS_PATH'), "db/sqlite/mysql2sqlite/mysql2sqlite")

#List of tables to update
tables_to_sync = ['ERF_IDs', 'CyberShake_Sites', 'Ruptures', 'CyberShake_Site_Ruptures', 'Rupture_Variations', 'Rup_Var_Seeds', 'CyberShake_Runs']

remote_conn = pymysql.connect(host=remote_db, user='cybershk_ro', passwd='CyberShake2007', db='CyberShake')
local_conn = sqlite3.connect(local_db)

remote_cur = remote_conn.cursor()
local_cur = local_conn.cursor()

local_tables_list = []
for r in local_cur.execute("select name from sqlite_master where type='table'").fetchall():
	local_tables_list.append(r[0])
print(local_tables_list)

for table in tables_to_sync:
	print("Checking table %s." % table)
	#Check number of entries - if it's different, drop the local version and recreate
	query = 'select count(*) from %s' % table
	remote_cur.execute(query)
	remote_count = int(remote_cur.fetchone()[0])
	if table in local_tables_list:
		local_cur.execute(query)
		local_count = int(local_cur.fetchone()[0])
	else:
		#This table doesn't exist locally, so set local_count to -1 to force an update
		local_count = -1
	update_table = False
	if remote_count!=local_count:
		print("Remote table %s has %d entries, local table %s has %d entries." % (table, remote_count, table, local_count))
		update_table = True
	elif force_update==True:
		print("Force was selected, will update table %s anyway." % table)
		update_table = True
	if update_table==True:
		print("Recreating local table.")
		
		#Drop local table
		if local_count>-1:
			drop_cmd = 'drop table %s' % table
			local_cur.execute(drop_cmd)
	
		#Get remote table with mysqldump
		mysqldump_cmd = 'mysqldump --single-transaction --host=%s --user=cybershk_ro --password=CyberShake2007 CyberShake %s > /tmp/%s_dump.sql' % (remote_db, table, table)
		rc = os.system(mysqldump_cmd)
		if rc!=0:
			print("Error running command %s, aborting." % mysqldump_cmd)
			sys.exit(rc)

		#Convert dump string into sqlite
		convert_cmd = '%s /tmp/%s_dump.sql > /tmp/%s_dump.sqlite' % (my2lite_path, table, table)
		rc = os.system(convert_cmd)
		if rc!=0:
			print("Error running command %s, aborting." % convert_cmd)
			sys.exit(rc)

		#Load new table
		with open('/tmp/%s_dump.sqlite' % table, 'r') as fp:
			table_load_string = fp.read()
			fp.close()
		local_cur.executescript(table_load_string)
		local_conn.commit()
remote_conn.close()
local_conn.close()

