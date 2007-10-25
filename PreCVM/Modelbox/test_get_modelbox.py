#! /usr/bin/python

"""
These are acceptance tests for the get_modelbox.py program.
"""
import unittest
import filecmp
import os

class Test_get_modelbox(unittest.TestCase):
  """ Acceptance Test for get_modelbox.py"""

  def test_pas(self):
    """Test new results against Site PAS"""

    ref_file = "./reference/PAS.modelbox"
    new_file = "./PAS.test_result"

    sstring = "./get_modelbox.py PAS "+new_file

    res = os.system(sstring)
    if (res/256) != 0:
      self.fail("Exit return error from get modelbox")

    self.failIf(filecmp.cmp(new_file,ref_file) == False,
	'new modelbox does not match existing modelbox.')
   
    rmstring = "rm " +new_file
    os.system(rmstring)

  def test_usc(self):
    """Test new results against Site USC"""

    ref_file = "./reference/USC.modelbox"
    new_file = "./USC.test_result"

    sstring = "./get_modelbox.py USC "+new_file

    res = os.system(sstring)
    if (res/256) != 0:
      self.fail("Exit return error from get modelbox")

    self.failIf(filecmp.cmp(new_file,ref_file) == False,
	'new modelbox does not match existing modelbox.')

    rmstring = "rm " +new_file
    os.system(rmstring)

  def test_ffi(self):
    """Test new results against Site FFI"""

    ref_file = "./reference/FFI.modelbox"
    new_file = "./FFI.test_result"

    sstring = "./get_modelbox.py FFI "+new_file

    res = os.system(sstring)
    if (res/256) != 0:
      self.fail("Exit return error from get modelbox")

    self.failIf(filecmp.cmp(new_file,ref_file) == False,
	'new modelbox does not match existing modelbox.')

    rmstring = "rm " +new_file
    os.system(rmstring)

  def test_smca(self):
    """Test new results against Site SMCA"""

    ref_file = "./reference/SMCA.modelbox"
    new_file = "./SMCA.test_result"

    sstring = "./get_modelbox.py SMCA "+new_file

    res = os.system(sstring)
    if (res/256) != 0:
      self.fail("Exit return error from get modelbox")

    self.failIf(filecmp.cmp(new_file,ref_file) == False,
	'new modelbox does not match existing modelbox.')

    rmstring = "rm " +new_file
    os.system(rmstring)

  def test_ccp(self):
    """Test new results against Site CCP"""

    ref_file = "./reference/CCP.modelbox"
    new_file = "./CCP.test_result"

    sstring = "./get_modelbox.py CCP "+new_file

    res = os.system(sstring)
    if (res/256) != 0:
      self.fail("Exit return error from get modelbox")

    self.failIf(filecmp.cmp(new_file,ref_file) == False,
	'new modelbox does not match existing modelbox.')

    rmstring = "rm " +new_file
    os.system(rmstring)

if __name__ == '__main__':
  unittest.main()
