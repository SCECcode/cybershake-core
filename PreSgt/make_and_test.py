#!/usr/bin/python

import sys
import subprocess
import os

os.chdir('src')
exitcode = subprocess.call('make clean')
if (exitcode <> 0):
	print "Make failed.\n"
	exit(1)

print "Make succeeded.\n"
os.chdir("..")
subprocess.call("./test_presgt.py")
