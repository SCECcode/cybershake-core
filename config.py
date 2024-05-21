#!/usr/bin/env python3

import os
import sys

vars = {}

def readCfg():
	try:
		filename ='%s/cybershake.cfg' % os.path.dirname(__file__)
		cfg = open(filename)
	except IOError:
		print("%s not found.\n" % filename)
		sys.exit(-2)
	cfgContents = cfg.readlines()
	for line in cfgContents:
		pieces = line.split('=')
		vars[pieces[0].strip()] = pieces[1].strip()

def getProperty(property):
	if len(vars)==0:
		readCfg()
	try:
		propertyVal = vars[property]
	except KeyError:
		raise KeyError("No %s found in cybershake.cfg.\n" % property)
	return propertyVal

def getJobID():
	jobid = os.getenv('PBS_JOBID')
	if jobid==None:
		jobid = os.getenv('JOB_ID')
	return jobid

	
