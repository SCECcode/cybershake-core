#!/usr/bin/env python

import MySQLdb
import sys
import os
import time

full_path = os.path.abspath(sys.argv[0])
path_add = os.path.dirname(os.path.dirname(full_path))

sys.path.append(path_add)

import config

def getSiteCoords(site):
	print "Retrieving coordinates from the database for %s.\n" % site
	#host = "focal.usc.edu"
	host = "moment.usc.edu"
	user = "cybershk_ro"
	passwd = "CyberShake2007"
	db = "CyberShake"
	try:
		cursor = MySQLdb.connect(host, user, passwd, db).cursor()
		sql_string = 'select CS_Site_Lat, CS_Site_Lon from CyberShake_Sites where CS_Short_Name="%s"' % site
		cursor.execute(sql_string)
		return cursor.fetchone()	
	except MySQLdb.OperationalError,e:
		print e
		sys.exit(-1)
	except MySQLdb.DatabaseError,e:
	        print e
	        sys.exit(-2)

def genFdloc(outputName, site, mlon, mlat, cordfileName):
    '''Duplicates gen_fdsrcloc.csh: generates an fdloc file, which contains the (X,Y) grid coordinates which are nearest to the site'''
    print "Generating %s.fdloc.\n" % site
    cordfile = open(cordfileName)
    minDistance = 1000.0
    minLoc = [-1, -1]
    i = 0
    for line in cordfile:
	    i = i+1
	    if i % 100000==0:
	        print '%d' % i
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


def genFaultList(outputName, site, erf_id):
    '''Copies the functionality of gen_faultlist.csh:  it serves as a wrapper to CreateFaultList.java, which queries the database to generate a list of ruptures which are applicable for the given site.'''
    print "Generating fault list for %s.\n" % site
    #command = 'java -classpath .:%s:%s/faultlist/mysql-connector-java-5.0.5-bin.jar faultlist/CreateFaultList %s %s %s %s' % (sys.path[0], sys.path[0], site, erf_id, PATH_TO_RUPTURE_VARIATIONS, outputName)
    command = '%s/faultlist_py/CreateFaultList.py %s %s %s %s' % (sys.path[0], site, erf_id, PATH_TO_RUPTURE_GEOMETRIES, outputName)
    print command
    returnCode = os.system(command)
    if returnCode!=0:
	sys.exit((returnCode >> 8) & 0xFF)

def genRadiusFile(radiusfile):
	'''Duplicates the functionality of part of gen_sgtgrid.csh:  it writes the adaptive mesh info to a radius file as part of generating the cordfile.'''

	print "Generating %s.\n" % radiusfile
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

	print "Generating %s.cordfile.\n" % site

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
                        print "Unable to read nodefile %s" % (node_file) 
                        sys.exit(1)
        elif (MPI_CMD == "aprun"):
		np = int(os.environ["PBS_NUM_NODES"])*4
		#No more than 32 cores)
		np = min(np, 32)
                MPI_CMD = "%s -n %d" % (MPI_CMD, np)

	command = '%s %s/bin/gen_sgtgrid nx=%d ny=%d nz=%d h=%f xsrc=%d ysrc=%d ixmin=%d ixmax=%d iymin=%d iymax=%d izstart=%d izmax=%d radiusfile=%s outfile=%s modellon=%f modellat=%f modelrot=%f faultlist=%s' % (MPI_CMD, sys.path[0], ns[0], ns[1], ns[2], HH, src[0], src[1], IX_MIN, IX_MAX, IY_MIN, IY_MAX, IZ_START, IZ_MAX, radiusfile, outputFile, mlon, mlat, mrot, faultlist)
	#cmdFile = open("command.txt", "w")
	#cmdFile.write(command)
	#cmdFile.flush()
	#cmdFile.close()
	print command
	startTime = time.time()
	returnCode = os.system(command)
	print "Elapsed time: %f\n" % (time.time()-startTime)
	if returnCode!=0:
		sys.exit((returnCode >> 8) & 0xFF)

RUPTURE_ROOT = config.getProperty('RUPTURE_ROOT')

if len(sys.argv) < 11:
    print 'Usage: ./presgt.py <site> <erf_id> <modelbox> <gridout> <model_coords> <fdloc> <faultlist> <radiusfile> <sgtcords> <spacing> [frequency]'
    print 'Example: ./presgt.py USC 33 USC.modelbox gridout_USC model_coords_GC_USC USC.fdloc USC.faultlist USC.radiusfile USC.cordfile 200.0 0.1'
    sys.exit(1)

site = sys.argv[1]
modelbox = sys.argv[3]
gridout = sys.argv[4]
cordfileName = sys.argv[5]

fdlocFileName = sys.argv[6]
faultlistFileName = sys.argv[7]
radiusFileName = sys.argv[8]
sgtcordFileName = sys.argv[9]

spacing = float(sys.argv[10])

frequency = 0.5
if len(sys.argv)==12:
	frequency = float(sys.argv[11])

erf_id = sys.argv[2]

PATH_TO_RUPTURE_GEOMETRIES = "%s/Ruptures_erf%s" % (RUPTURE_ROOT, erf_id)

input = open(modelbox)
modelboxContents = input.readlines()
input.close()
site = modelboxContents[0].strip()

siteCoords = getSiteCoords(site)
print siteCoords
siteLat = float(siteCoords[0])
siteLon = float(siteCoords[1])

modelTokens = modelboxContents[4].split()
mlon = (float)(modelTokens[1])
mlat = (float)(modelTokens[3])
mrot = (float)(modelTokens[5])

src = genFdloc(fdlocFileName, site, siteLon, siteLat, cordfileName)
genFaultList(faultlistFileName, site, erf_id)

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

