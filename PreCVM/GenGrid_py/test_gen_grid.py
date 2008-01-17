#!/usr/bin/python

import unittest
import os
import sys
import gen_grid
import compare

class Test_get_grid(unittest.TestCase):
    """ Acceptance Test for gen_grid.py"""	
	
    def testPAS(self):
        referenceFiles = ["reference/PAS/gridout_PAS", "reference/PAS/gridfile_PAS", "reference/PAS/model_bounds_GC_PAS", "reference/PAS/model_coords_GC_PAS", "reference/PAS/model_params_GC_PAS"]

        testFiles = ["gridout_PAS", "gridfile_PAS", "model_bounds_GC_PAS", "model_coords_GC_PAS", "model_params_GC_PAS"]
		
		
        gen_grid.genGrid("reference/PAS/PAS.modelbox", testFiles[1], testFiles[0], testFiles[3], testFiles[4], testFiles[2])
        for i in range(5):
            self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "File " + testFiles[i] + " doesn't match " + referenceFiles[i] + ".")
        
        for file in testFiles:
            rmstring = "rm " + file
            os.system(rmstring)
        
    
    def testFFI(self):
        referenceFiles = ["reference/FFI/gridout_FFI", "reference/FFI/gridfile_FFI", "reference/FFI/model_bounds_GC_FFI", "reference/FFI/model_coords_GC_FFI", "reference/FFI/model_params_GC_FFI"]

        testFiles = ["gridout_FFI", "gridfile_FFI", "model_bounds_GC_FFI", "model_coords_GC_FFI", "model_params_GC_FFI"]
        
        gen_grid.genGrid("reference/FFI/FFI.modelbox", testFiles[1], testFiles[0], testFiles[3], testFiles[4], testFiles[2])
        for i in range(5):
            self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "File " + testFiles[i] + " doesn't match " + referenceFiles[i] + ".")
        
        for file in testFiles:
            rmstring = "rm " + file
            os.system(rmstring)
        
    
    def testCCP(self):
        referenceFiles = ["reference/CCP/gridout_CCP", "reference/CCP/gridfile_CCP", "reference/CCP/model_bounds_GC_CCP", "reference/CCP/model_coords_GC_CCP", "reference/CCP/model_params_GC_CCP"]

        testFiles = ["gridout_CCP", "gridfile_CCP", "model_bounds_GC_CCP", "model_coords_GC_CCP", "model_params_GC_CCP"]
        
        gen_grid.genGrid("reference/CCP/CCP.modelbox", testFiles[1], testFiles[0], testFiles[3], testFiles[4], testFiles[2])
        for i in range(5):
            self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "File " + testFiles[i] + " doesn't match " + referenceFiles[i] + ".")
        
        for file in testFiles:
            rmstring = "rm " + file
            os.system(rmstring)
        
    
    def testSMCA(self):
        referenceFiles = ["reference/SMCA/gridout_SMCA", "reference/SMCA/gridfile_SMCA", "reference/SMCA/model_bounds_GC_SMCA", "reference/SMCA/model_coords_GC_SMCA", "reference/SMCA/model_params_GC_SMCA"]

        testFiles = ["gridout_SMCA", "gridfile_SMCA", "model_bounds_GC_SMCA", "model_coords_GC_SMCA", "model_params_GC_SMCA"]
        
        gen_grid.genGrid("reference/SMCA/SMCA.modelbox", testFiles[1], testFiles[0], testFiles[3], testFiles[4], testFiles[2])
        for i in range(5):
            self.failIf(not compare.compare(referenceFiles[i], testFiles[i], 0.00011), "File " + testFiles[i] + " doesn't match " + referenceFiles[i] + ".")
        
        for file in testFiles:
            rmstring = "rm " + file
            os.system(rmstring)
    

if __name__ == '__main__':
  unittest.main()
