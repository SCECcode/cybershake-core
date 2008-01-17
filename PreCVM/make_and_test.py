#!/usr/bin/python

import os
import sys


def makeGetpar():
	os.chdir("../Getpar/getpar/src")
	os.system("make clean")
	exitcode = os.system("make")
	if (exitcode <> 0):
		print "Make getpar failed.\n"
		sys.exit(1)
	else:
		print "Make getpar succeeded.\n"
		os.chdir("../../../PreCVM")
	

def makeModelbox():
	os.chdir("Modelbox/src")
	os.system("make clean")
	exitcode = os.system("make")
	if (exitcode <> 0):
		print "Make modelbox failed.\n"
		sys.exit(2)
	else:
		print "Make modelbox succeeded.\n"
	os.chdir("../..")

def makeGenGrid():
	os.chdir(sys.path[0]+"/GenGrid_py/src")
	os.system("make clean")
	exitcode = os.system("make")
	if (exitcode <> 0):
		print "Make gen_grid failed.\n"
		sys.exit(3)
	else:
		print "Make gen_grid succeeded, testing.\n"
	os.chdir("../..")


def runTest():
	os.system("./test_pre_cvm.py");


makeGetpar()
makeModelbox()
makeGenGrid()
runTest()
