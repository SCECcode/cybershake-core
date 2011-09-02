#!/usr/bin/env python

import os
import sys

vars = {}

def readCfg():
	try:
		cfg_file = os.environ['CYBERSHAKE_CONFIG']
	except KeyError:
		#default to cybershake.cfg in same directory
		cfg_file = '%s/cybershake.cfg' % os.path.dirname(__file__)
	try:
		cfg = open(cfg_file)
        except IOError:
                print "%s not found.\n" % cfg_file
		sys.exit(-2)
        cfgContents = cfg.readlines()
	for line in cfgContents:
        	pieces = line.split('=')
		if pieces[0][0]=='#':
			continue
        	vars[pieces[0].strip()] = pieces[1].strip()

def getProperty(property):
	if len(vars)==0:
		readCfg()
	try:
		propertyVal = vars[property]
	except KeyError:
		print "No %s found in cybershake.cfg.\n" % property
		sys.exit(-1)
	return propertyVal

def getJobID():
	jobid = os.getenv('PBS_JOBID')
	if jobid==None:
		jobid = os.getenv('JOB_ID')
	return jobid

	
