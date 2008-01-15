#!/usr/bin/python

import unittest
import subprocess
import sys
import os
import compare

class Test_presgt(unittest.TestCase):
	'''Acceptance tests for PreSGT steps.'''
	ERF_ID=33

	def testSites(self):
		siteList = ['USC', 'PAS', 'SMCA', 'CCP', 'FFI']
		for site in siteList:
			referenceFiles = ['reference/FdLocs/%s.fdloc' % site, 'reference/RadiusFile/%s.radiusfile' % site, 'reference/FaultList/%s.faultlist' % site, 'reference/SgtCords/%s.cordfile' % site]
			testFiles = ['%s.fdloc' % site, '%s.radiusfile' % site, '%s.faultlist' % site, '%s.cordfile' % site]
			subprocess.call("./presgt.py %s %d" % (site, ERF_ID))
			for i in range(5):
				self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "File " + testFiles[i] + " doesn't match " + referenceFiles[i] + ".")

			for file in testFiles:
				os.system("rm %s" % file)
			

if __name__ == '__main__':
	unittest.main()
