#!/usr/bin/python

import sys
import subprocess
import os

os.chdir('src')
subprocess.call(['make', 'clean'])
exitcode = subprocess.call('make')
if (exitcode <> 0):
	print "Make failed.\n"
	exit(1)

print "Make succeeded.\n"
os.chdir("..")
subprocess.call("./test_presgt.py")
