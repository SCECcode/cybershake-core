#!/bin/csh

if ( $# != 7 ) then
	echo "Usage: `basename $0` SITE NX NY NZ MODEL_LAT MODEL_LON MODEL_ROT"
	echo ""
	echo "SITE: the site id (e.g. LADT)"
	echo "NX:"
	echo "NY:"
	echo "NZ:"
	echo "MODEL_LAT:"
	echo "MODEL_LON:"
	echo "MODEL_ROT:"
endif


set SITE = $1
set NX = $2
set NY = $3
set NZ = $4
set MODEL_LAT = $5
set MODEL_LON = $6
set MODEL_ROT = $7

set HH = 0.2

set IX_MIN = 20
set IX_MAX = `echo $NX | gawk '{printf "%d\n",$1-20;}'`

set IY_MIN = 20
set IY_MAX = `echo $NY | gawk '{printf "%d\n",$1-20;}'`

set IZ_START = 1
set IZ_MAX = 180

set RLEV = ( 10.0 50.0 100.0 1000.0 )
set RINC = (   10   15    25     50 )

set ZLEV = (  0.2 5.0 24.0 60.0 )
set ZINC = (    1   5   10   25 )

set RADIUS_FILE = ../../data/SgtInfo/RadiusFile/${SITE}.radiusfile

echo $#RLEV > $RADIUS_FILE
echo $RLEV >> $RADIUS_FILE
echo $RINC >> $RADIUS_FILE
echo $#ZLEV >> $RADIUS_FILE
echo $ZLEV >> $RADIUS_FILE
echo $ZINC >> $RADIUS_FILE

chgrp lc_pjm $RADIUS_FILE
chmod g+rw $RADIUS_FILE

set FAULTLIST = ../../data/SgtInfo/FaultList/${SITE}.faultlist

set SRC_CORDS = ../../data/SgtInfo/FdLocs/${SITE}.fdloc
set XSRC = `gawk '{printf "%d\n",$1;}' $SRC_CORDS `
set YSRC = `gawk '{printf "%d\n",$2;}' $SRC_CORDS `
echo $XSRC $YSRC

set OUTDIR = ../../data/SgtInfo/SgtCords/
mkdir -p $OUTDIR
set OUTFILE = $OUTDIR/${SITE}.cordfile
rm -f $OUTFILE

bin/gen_sgtgrid nx=$NX ny=$NY nz=$NZ h=$HH xsrc=$XSRC ysrc=$YSRC \
            ixmin=$IX_MIN ixmax=$IX_MAX iymin=$IY_MIN iymax=$IY_MAX \
	    izstart=$IZ_START izmax=$IZ_MAX \
	    radiusfile=$RADIUS_FILE outfile=$OUTFILE \
	    modellon=$MODEL_LON modellat=$MODEL_LAT modelrot=$MODEL_ROT \
	    faultlist=$FAULTLIST

set RC = $?

chgrp lc_pjm $OUTFILE
chmod g+w $OUTFILE

if ( $RC != 0 ) exit $RC
