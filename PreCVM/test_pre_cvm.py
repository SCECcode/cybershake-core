#!/usr/bin/python

import unittest
import subprocess
import filecmp
import compare

class Test_pre_cvm(unittest.TestCase):
	'''Acceptance test for pre_cvm.py'''

	def testAll(self):
		siteList=['PAS', 'SMCA', 'CCP', 'FFI']
		for site in siteList:
			modelboxRef = 'Modelbox/reference/%s.modelbox' % site
			gridfileRef = 'GenGrid_py/reference/%s/gridfile_%s' % (site, site)
			gridoutRef = 'GenGrid_py/reference/%s/gridout_%s' % (site, site)
			coordsfileRef = 'GenGrid_py/reference/%s/model_coords_GC_%s' % (site, site)
			paramsfileRef = 'GenGrid_py/reference/%s/model_params_GC_%s' % (site, site)
			boundsfileRef = 'GenGrid_py/reference/%s/model_bounds_GC_%s' % (site, site)
			referenceFiles = [gridfileRef, gridoutRef, coordsfileRef, paramsfileRef, boundsfileRef]

			modelbox = '%s.modelbox' % site
			gridfile = 'gridfile_%s' % site
			gridout = 'gridout_%s' % site
			coordsfile = 'model_coords_%s' % site
			paramsfile = 'model_params_%s' % site
			boundsfile = 'model_bounds_%s' % site
			testFiles = [gridfile, gridout, coordsfile, paramsfile, boundsfile]
			
			subprocess.call(['./pre_cvm.py', site, modelbox, testFiles[0], testFiles[1], testFiles[2], testFiles[3], testFiles[4]])

			self.failIf(not filecmp.cmp(modelboxRef, modelbox), "%s doesn't match %s." % (modelboxRef, modelbox))
			subprocess.call(["rm",modelbox])
			
			for i in range(len(referenceFiles)):
				self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "%s doesn't match %s." % (referenceFiles[i], testFiles[i]))

			for file in testFiles:
				subprocess.call(["rm",file])


if __name__=='__main__':
	unittest.main()
