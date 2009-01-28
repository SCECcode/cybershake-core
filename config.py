#!/usr/bin/env python

import os
import sys

vars = {}

def readCfg():
	cfg = open('cybershake.cfg')
        if cfg==None:
                print "%s not found." % os.path.join(os.getcwd(), 'cybershake.cfg')
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
		print "No %s found in cybershake.cfg.\n" % property
		sys.exit(-1)
	return propertyVal

