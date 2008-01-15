#!/usr/bin/python

import unittest
import subprocess
import sys
import os
import compare

class Test_presgt(unittest.TestCase):
	'''Acceptance tests for PreSGT steps.'''

	def testSites(self):
		ERF_ID=33
		siteList = ['USC', 'PAS', 'SMCA', 'CCP', 'FFI']
		for site in siteList:
			referenceFiles = ['reference/SgtInfo/FdLocs/%s.fdloc' % site, 'reference/SgtInfo/RadiusFile/%s.radiusfile' % site, 'reference/SgtInfo/FaultList/%s.faultlist' % site, 'reference/SgtInfo/SgtCords/%s.cordfile' % site]
			testFiles = ['%s.fdloc' % site, '%s.radiusfile' % site, '%s.faultlist' % site, '%s.cordfile' % site]
			#args = "%s %d %s %s %s" % (site, ERF_ID, 'reference/ModelCords/' + site + '/' + site + '_modelbox', 'reference/ModelCords/' + site + '/gridout_' + site, 'reference/ModelCords/' + site + '/model_coords_GC_' + site)
			subprocess.call(["./presgt.py", site, "%d" % ERF_ID, 'reference/ModelParams/%s/%s.modelbox' % (site, site), 'reference/ModelParams/%s/gridout_%s' % (site, site), 'reference/ModelParams/%s/model_coords_GC_%s' % (site, site)])
			for i in range(5):
				self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "File " + testFiles[i] + " doesn't match " + referenceFiles[i] + ".")

			for file in testFiles:
				os.system("rm %s" % file)
			

if __name__ == '__main__':
	unittest.main()
