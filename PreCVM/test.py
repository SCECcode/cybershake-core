#!/usr/bin/python

import os
import sys


def doModelbox():
	os.chdir(sys.path[0]+"/Modelbox")
	exitcode = os.system("./test_get_modelbox.py")
	if (exitcode <> 0):
		print "Test modelbox failed.\n"
		return
	else:
		print "Modelbox testing successful!\n"

def doGenGrid():
	os.chdir(sys.path[0]+"/GenGrid_py")
	exitcode = os.system("./test_gen_grid.py")
	if (exitcode <> 0):
		print "Test genGrid failed.\n"
		return
	else:
		print "GenGrid testing successful!\n"


doModelbox()
doGenGrid()

