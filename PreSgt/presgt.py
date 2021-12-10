#!/usr/bin/env python

import pymysql
import sys
import os
import time
import sqlite3
import optparse

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

def getSiteCoords(site, host):
	print("Retrieving coordinates from the database for %s.\n" % site)
	user = "cybershk_ro"
	passwd = "CyberShake2007"
	db = "CyberShake"
	if host[0:9]=="sqlite://":
		cursor = sqlite3.connect(host[9:]).cursor()
	else:
		cursor = pymysql.connect(host, user, passwd, db).cursor()
	sql_string = 'select CS_Site_Lat, CS_Site_Lon from CyberShake_Sites where CS_Short_Name="%s"' % site
	cursor.execute(sql_string)
	return cursor.fetchone()	

def genFdloc(outputName, site, mlon, mlat, cordfileName):
    '''Duplicates gen_fdsrcloc.csh: generates an fdloc file, which contains the (X,Y) grid coordinates which are nearest to the site'''
    print("Generating %s.fdloc.\n" % site)
    cordfile = open(cordfileName)
    minDistance = 1000.0
    minLoc = [-1, -1]
    i = 0
    for line in cordfile:
	    i = i+1
	    if i % 100000==0:
	        print('%d' % i)
	    pieces = line.split()
	    distance = (float(pieces[0])-mlon)*(float(pieces[0])-mlon)+(float(pieces[1])-mlat)*(float(pieces[1])-mlat)
	    if distance < minDistance:
	       	minDistance = distance
	       	minLoc[0] = pieces[2]
	       	minLoc[1] = pieces[3]
    cordfile.close()
    fdlocFile = open(outputName, 'w')
    fdlocFile.write('%s %s\n' % (minLoc[0], minLoc[1]))
    fdlocFile.flush()
    fdlocFile.close()
    return [int(minLoc[0]), int(minLoc[1])]


def genFaultList(outputName, site, erf_id, rsqsim):
    '''Copies the functionality of gen_faultlist.csh:  it serves as a wrapper to CreateFaultList.java, which queries the database to generate a list of ruptures which are applicable for the given site.'''
    print("Generating fault list for %s.\n" % site)
    #command = 'java -classpath .:%s:%s/faultlist/mysql-connector-java-5.0.5-bin.jar faultlist/CreateFaultList %s %s %s %s' % (sys.path[0], sys.path[0], site, erf_id, PATH_TO_RUPTURE_VARIATIONS, outputName)
    rsqsim_str = ""
    if rsqsim:
        rsqsim_str = "--rsqsim"
    command = '%s/faultlist_py/CreateFaultList.py --site %s --erf_id %s --rup_path %s --output %s --server %s %s' % (sys.path[0], site, erf_id, PATH_TO_RUPTURE_GEOMETRIES, outputName, server, rsqsim_str)
    print(command)
    returnCode = os.system(command)
    if returnCode!=0:
        sys.exit((returnCode >> 8) & 0xFF)

def genRadiusFile(radiusfile):
	'''Duplicates the functionality of part of gen_sgtgrid.csh:  it writes the adaptive mesh info to a radius file as part of generating the cordfile.'''

	print("Generating %s.\n" % radiusfile)
	#magic constants for the adaptive mesh
	RLEV = [10.0, 50.0, 100.0, 1000.0]
	RINC = [10, 15, 25, 50]
	ZLEV = [0.2, 5.0, 24.0, 60.0]
	ZINC = [1, 5, 10, 25]

	output = open(radiusfile, 'w')
	output.write('%d\n' % len(RLEV))
	for r in RLEV:
		output.write('%f ' % r)
	output.write('\n')
	for r in RINC:
		output.write('%d ' % r)
	output.write('\n')
	output.write('%d\n'% len(ZLEV))
	for z in ZLEV:
		output.write('%f ' % z)
	output.write('\n')
	for z in ZINC:
		output.write('%d ' % z)
	output.write('\n')
	output.flush()
	output.close()


def genSgtGrid(outputFile, site, ns, src, mlon, mlat, mrot, faultlist, radiusfile, frequency, spacing):
	'''Copies the functionality of gen_sgtgrid.csh:  it generates a list of grid points that the SGTs should be saved for.  This includes an adaptive mesh as well as the locations of the ruptures.  Essentially this serves as a wrapper for gen_sgtgrid.c.'''
	#magic constants
	HH = 0.1/frequency
	if spacing>0:
		HH = spacing
	IX_MIN = 20
	IX_MAX = ns[0]-20
	IY_MIN = 20
	IY_MAX = ns[1]-20
	IZ_START = 1
	IZ_MAX = ns[2]-20

	#outfile = '../../data/SgtInfo/SgtCords/%s.cordfile' % site
	#faultlist = '../../data/SgtInfo/FaultList/%s.faultlist' % site

	print("Generating %s.cordfile.\n" % site)

	#split into subfiles
	MPI_CMD = config.getProperty('MPI_CMD')

	if (MPI_CMD == "mpirun"):
                try:
                        node_file = os.environ["PBS_NODEFILE"]
                        num_nodes = 0
                        f = open(node_file)
                        lines = f.readlines()
                        num_nodes = len(lines)
                        MPI_CMD = "%s -np %d -machinefile %s" % \
			    (MPI_CMD, num_nodes, node_file)
                except:
                        print("Unable to read nodefile %s" % (node_file))
                        sys.exit(1)
	elif (MPI_CMD == "aprun"):
		np = int(os.environ["PBS_NUM_NODES"])*4
		#No more than 32 cores)
		np = min(np, 32)
		MPI_CMD = "%s -n %d -N 4" % (MPI_CMD, np)
	elif (MPI_CMD == "jsrun"):
		num_res_sets = int(os.environ['LSB_DJOB_NUMPROC'])-1
		#No more than 32 cores
		num_res_sets = min(num_res_sets, 32)
		MPI_CMD = "%s -a 1 -c 1 -r %d -n %d" % (MPI_CMD, num_res_sets, num_res_sets)
	command = '%s %s/bin/gen_sgtgrid nx=%d ny=%d nz=%d h=%f xsrc=%d ysrc=%d ixmin=%d ixmax=%d iymin=%d iymax=%d izstart=%d izmax=%d radiusfile=%s outfile=%s modellon=%f modellat=%f modelrot=%f faultlist=%s' % (MPI_CMD, sys.path[0], ns[0], ns[1], ns[2], HH, src[0], src[1], IX_MIN, IX_MAX, IY_MIN, IY_MAX, IZ_START, IZ_MAX, radiusfile, outputFile, mlon, mlat, mrot, faultlist)
	#cmdFile = open("command.txt", "w")
	#cmdFile.write(command)
	#cmdFile.flush()
	#cmdFile.close()
	print(command)
	startTime = time.time()
	returnCode = os.system(command)
	print("Elapsed time: %f\n" % (time.time()-startTime))
	if returnCode!=0:
		sys.exit((returnCode >> 8) & 0xFF)

RUPTURE_ROOT = config.getProperty('RUPTURE_ROOT')

parser = optparse.OptionParser()
parser.add_option("--site", dest="site", action="store", help="Site name")
parser.add_option("--erf_id", dest="erf_id", action="store", type="int", help="ERF ID")
parser.add_option("--modelbox", dest="modelbox", action="store", help="Path to modelbox file (input)")
parser.add_option("--gridout", dest="gridout", action="store", help="Path to gridout file (input)")
parser.add_option("--coordfile", dest="coordsfile", action="store", help="Path to model_coords file (input)")
parser.add_option("--fdloc", dest="fdlocfile", action="store", help="Path to fdloc file (output)")
parser.add_option("--faultlist", dest="faultlistfile", action="store", help="Path to faultlist file (otuput)")
parser.add_option("--radiusfile", dest="radiusfile", action="store", help="Path to radiusfile (output)")
parser.add_option("--sgtcords", dest="sgtcordsfile", action="store", help="Path to SGT coords file (output)")
parser.add_option("--spacing", dest="spacing", action="store", type="float", help="Mesh spacing, in km")
parser.add_option("--frequency", dest="frequency", action="store", default=0.5, type="float", help="Override default frequency of 0.5 Hz")
parser.add_option("--rsqsim", dest="rsqsim", action="store_true", default=False, help="Assumes RSQSim-formatted rupture geometry files.  Default is false.")
parser.add_option("--server", dest="server", action="store", default="moment.usc.edu", help="Path to server (default is moment.usc.edu).  Can specify SQLite file with sqlite://<file>")

(option, args) = parser.parse_args()

site = option.site
erf_id = option.erf_id
modelbox = option.modelbox
gridout = option.gridout
cordfileName = option.coordsfile

fdlocFileName = option.fdlocfile
faultlistFileName = option.faultlistfile
radiusFileName = option.radiusfile
sgtcordFileName = option.sgtcordsfile

spacing = option.spacing

if (site==None or erf_id==None or modelbox==None or gridout==None or cordfileName==None or fdlocFileName==None or faultlistFileName==None or radiusFileName==None or sgtcordFileName==None or spacing==None):
	print("All of site, erf_id, modelbox, gridout, coordfile, floc, faultlist, radiusfile, sgtcords, and spacing must be specified.")
	parser.print_help()
	sys.exit(1)

rsqsim = option.rsqsim
frequency = option.frequency

server = option.server
	
PATH_TO_RUPTURE_GEOMETRIES = "%s/Ruptures_erf%s" % (RUPTURE_ROOT, erf_id)

input = open(modelbox)
modelboxContents = input.readlines()
input.close()
site = modelboxContents[0].strip()

siteCoords = getSiteCoords(site, server)
siteLat = float(siteCoords[0])
siteLon = float(siteCoords[1])

modelTokens = modelboxContents[4].split()
mlon = (float)(modelTokens[1])
mlat = (float)(modelTokens[3])
mrot = (float)(modelTokens[5])

src = genFdloc(fdlocFileName, site, siteLon, siteLat, cordfileName)
genFaultList(faultlistFileName, site, erf_id, rsqsim)

genRadiusFile(radiusFileName)

input = open(gridout)
gridoutContents = input.readlines();
input.close()
ns = []
ns.append(int((gridoutContents[1].split("="))[1]))
ns.append(int((gridoutContents[1+ns[0]+2].split("="))[1]))
ns.append(int((gridoutContents[1+ns[0]+2+ns[1]+2].split("="))[1]))

genSgtGrid(sgtcordFileName, site, ns, src, mlon, mlat, mrot, faultlistFileName, radiusFileName, frequency, spacing)
sys.exit(0)

