#! /usr/bin/env python

"""
This program inputs the ruptures that are relevant to a particular 
site and outputs a fault centroid and associated
boundaries for the simulation box.
"""

#
# Read standard python startup files. Performed before local imports
# because site startup file may drop imported modules from the env
#
import os
filename = os.environ.get('PYTHONSTARTUP')
if filename and os.path.isfile(filename):
  execfile(filename)

import os
import sys
import string
import math
import MySQLdb

class faultinfo:
  """
    class to store fault information returned from database
  """
  def __init__(self):
    mnlat = 0.0
    mnlon = 0.0
    mxlat = 0.0
    mxlon = 0.0

class boundaryinfo:
  """
    class to store fault information returned from database
  """
  def __init__(self):
    lat = 0.0
    lon = 0.0
    min = 0.0
    max = 0.0

def opendb(h,pt,us,pwd,d):
  try:
    connection = MySQLdb.connect(host=h,port=pt,user=us,passwd=pwd,db=d)
    return connection
  except MySQLdb.OperationalError,message:
    errorMessage = "Error %d:\n%s" % (message [ 0],message[1])
    print errorMessage
    sys.exit(-1)

def closedb(conn):
  try:
    conn.close()
    return	
  except MySQLdb.OperationalError,message:
    errorMessage = "Error %d:\n%s" % (message [ 0],message[1])
    print errorMessage
    sys.exit(-1)

#
# Special reserved magic numbers used in this script
# Also database password and 
#
pi= 3.14159
kplat = 111.19
model_rot = -55.0
bound_pad = 30.0
xmin = 1.0e15
xmax = -1.0e15
ymin = 1.0e15
ymax = -1.0e15
host = "focal.usc.edu"
#username = "scottcal"
#host = "hpc-scec.usc.edu"
username = "cybershk_ro"
pwd = "CyberShake2007"
port = 3306
db="CyberShake"

if len(sys.argv) < 4:
  print "Syntax: get_modelbox.py SITE_Name ERF_ID Outfile_name"
  print "Example: get_modelbox.py USC 34 ./usc_modelbox.txt"
  sys.exit()

site = sys.argv[1]
erf = sys.argv[2]
outfile = sys.argv[3]

gpu = False
if len(sys.argv) >= 5 and sys.argv[4]=="gpu":
	print "GPU mode."
	gpu = True

f = open(outfile,"w")

f.write("%s\n"%(site))

# Open CyberShake DB
conn = opendb(host,port,username,pwd,db)
cur = conn.cursor()

sql_string = "select * from CyberShake_Sites"
try:
  cur.execute(sql_string)
except MySQLdb.DatabaseError,e:
  print e
  sys.exit(-1)
  
station_set = cur.fetchall()

found = 0
for stas in station_set:
  if site == stas[2]:
    siteid=stas[0]
    found = 1
    break

if found !=1:
   print "Unable to find station %s in DB. Exiting\n" % site
   sys.exit(-1)
 
sql_string = "select distinct Ruptures.Source_ID,Ruptures.Rupture_ID,Source_Name,Mag,Prob,Start_Lat,Start_Lon,End_Lat,End_Lon from CyberShake_Sites,Ruptures,CyberShake_Site_Ruptures where Ruptures.ERF_ID=%s and Ruptures.ERF_ID=CyberShake_Site_Ruptures.ERF_ID and CyberShake_Site_Ruptures.CS_Site_ID=%s and CyberShake_Site_Ruptures.Source_ID=Ruptures.Source_ID and CyberShake_Site_Ruptures.Rupture_ID=Ruptures.Rupture_ID order by Mag desc"%(erf,stas[0])
try:
  cur.execute(sql_string)
except MySQLdb.DatabaseError,e:
  print e
  sys.exit(-1)

count = cur.rowcount

rup_set = cur.fetchall()

#
# Close db, connections, and results sets
#
cur.close()
closedb(conn)

#
# Initialize Fault list
#
faults = []

#
# Transfer result set to faults list
#

for all_rups in rup_set:
  finfo = faultinfo()
  finfo.mnlon = all_rups[6]
  finfo.mnlat = all_rups[5]
  finfo.mxlon = all_rups[8]
  finfo.mxlat = all_rups[7]
  faults.append(finfo)

#
# Initialize with Static parameters defined at top of file
#
mrot = model_rot
bpad = bound_pad

#
# Initialize the cumulative counters, and set lon,lat lists to empty
#
clon = 0.0
clat = 0.0
np = 0
lon = []
lat = []

for flts in faults:
  np = np+1
  lon.append(flts.mnlon)
  lat.append(flts.mnlat)
  clon = clon+flts.mnlon
  clat = clat+flts.mnlat

  np = np+1
  lon.append(flts.mxlon)
  lat.append(flts.mxlat)
  clon = clon+flts.mxlon
  clat = clat+flts.mxlat
  
clon = clon/(1.0*np)
clat = clat/(1.0*np)

kplon = kplat*math.cos((pi*clat)/180.0)
cosR = math.cos(pi*(90.0+mrot)/180.0)
sinR = math.sin(pi*(90.0+mrot)/180.0)

xm = []
ym = []

for i in range(0,np-1):
  n = (lat[i]-clat)*kplat
  e = (lon[i]- clon)*kplon
  xm.append((e*sinR) + (n*cosR))
  ym.append((e*cosR) - (n*sinR))
  if xm[i]<xmin:
    xmin=xm[i]
  if xm[i]>xmax:
    xmax=xm[i]
    
  if ym[i]<ymin:
    ymin=ym[i]
  if ym[i]>ymax:
    ymax=ym[i]

xmax = xmax+bpad
xmin = xmin-bpad
ymax = ymax+bpad
ymin = ymin-bpad

xlen=xmax-xmin
ylen=ymax-ymin

x0=xmax-0.5*xlen
y0=ymax-0.5*ylen

n= (-y0*sinR) + (x0*cosR)
e = (y0*cosR) + (x0*sinR)

mlon=clon + (e/kplon)
mlat=clat + (n/kplat)

#Round up to nearest 10km
xlrnd=10*int((xlen/10.0) + 0.5)
ylrnd=10*int((ylen/10.0) + 0.5)

if gpu:
        #Round X to the nearest multiple of 8 km, Y to 4 km, so that the # of grid points will be a multiple of 80 in X and 40 in Y (for 1 Hz)
        print "Old lengths: %d, %d" % (xlrnd, ylrnd)
        xlrnd += (-1*xlrnd) % 8
        ylrnd += (-1*ylrnd) % 4
        print "New lengths: %d, %d" % (xlrnd, ylrnd)


f.write("APPROXIMATE CENTROID:\n")
f.write(" clon=%12.5f clat=%12.5f\n"%(clon,clat))
f.write("MODEL PARAMETERS:\n")
f.write(" mlon= %12.5f mlat= %12.5f mrot= %6.1f xlen= %12.5f ylen= %12.5f\n"%(mlon,mlat,mrot,xlrnd,ylrnd))
f.write("MODEL CORNERS:\n")
systring ="echo 0.0 0.0 %f 0.0 %f %f 0.0 %f | bin/gcproj ref_lon=%9.5f ref_lat=%7.5f alpha=%f xlen=%f ylen=%f\n"%(xlrnd,xlrnd,ylrnd,ylrnd,mlon,mlat,mrot,xlrnd,ylrnd)
gcout = os.popen(systring)

vals = []

for x,lin in enumerate(gcout):
  val = lin.split()
  if x%2==0:
    mb = boundaryinfo()
    mb.lon=float(val[0])
    mb.lat=float(val[1])
  else:
    mb.min=float(val[0])
    mb.max=float(val[1])
    vals.append(mb)
  
if len(vals) != 4:
  print "gcproj did not return expected values"
  sys.exit(-1)

for vs in vals: 
  f.write("  %12.5f %12.5f (x= %10.3f y= %10.3f)\n"%(vs.lon,vs.lat,vs.min,vs.max))

f.close()
sys.exit(0)
