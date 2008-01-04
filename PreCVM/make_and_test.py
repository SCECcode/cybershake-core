#!/usr/bin/python

import os
import sys


def doGetpar():
	os.chdir("../Getpar/getpar/src")
	os.system("make clean")
	exitcode = os.system("make")
	if (exitcode <> 0):
		print "Make getpar failed.\n"
		exit(1)
	else:
		print "Make getpar succeeded.\n"
		os.chdir("../../../PreCVM")
	

def doModelbox():
	os.chdir("Modelbox/src")
	os.system("make clean")
	exitcode = os.system("make")
	if (exitcode <> 0):
		print "Make modelbox failed.\n"
		return
	else:
		print "Make modelbox succeeded, testing.\n"
	os.chdir("..")
	exitcode = os.system("./test_get_modelbox.py")
	if (exitcode <> 0):
		print "Test modelbox failed.\n"
		return
	else:
		print "Modelbox testing successful!\n"

def doGenGrid():
	os.chdir(sys.path[0]+"/GenGrid_py/src")
	os.system("make clean")
	exitcode = os.system("make")
	if (exitcode <> 0):
		print "Make gen_grid failed.\n"
		return
	else:
		print "Make gen_grid succeeded, testing.\n"
	os.chdir("..")
	exitcode = os.system("./test_gen_grid.py")
	if (exitcode <> 0):
		print "Test genGrid failed.\n"
		return
	else:
		print "GenGrid testing successful!\n"


doGetpar()
doModelbox()
doGenGrid()

