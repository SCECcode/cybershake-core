#!/bin/csh

set SITE = LADT

set OUTDIR = ../../data/SgtInfo/SgtCords/

mkdir -p $OUTDIR

set STATS = ( ${SITE} )

set HH = 0.2

set NX = 2100
set NY = 2850
set NZ = 200
#set SUFX = _sml-h${HH}

#set NX = 2500
#set NY = 2500
#set NZ = 250
#set SUFX = -h${HH}

set IX_MIN = 20
set IX_MAX = `echo $NX | gawk '{printf "%d\n",$1-20;}'`

set IY_MIN = 20
set IY_MAX = `echo $NY | gawk '{printf "%d\n",$1-20;}'`

set IZ_START = 1
set IZ_MAX = 180

set MODEL_LAT = 34.52629
set MODEL_LON = -118.
set MODEL_ROT = -90.0

set RLEV = ( 10.0 50.0 100.0 1000.0 )
set RINC = (   10   15    25     50 )

set ZLEV = (  5.0 24.0 60.0 )
set ZINC = (    5   10   25 )

#set RADIUS_FILE = sgt.radiusfile
set RADIUS_FILE = ../../data/SgtInfo/RadiusFile/${SITE}.radiusfile

echo $#RLEV > $RADIUS_FILE
echo $RLEV >> $RADIUS_FILE
echo $RINC >> $RADIUS_FILE
echo $#ZLEV >> $RADIUS_FILE
echo $ZLEV >> $RADIUS_FILE
echo $ZINC >> $RADIUS_FILE

set FAULTLIST = ../../data/SgtInfo/FaultList/${SITE}.faultlist

#set FAULTLIST = fault_list.file
#set RUPDIR = ../RupModel/RupturesFromUSC/RuptureVariations
#set RUPS = `\ls $RUPDIR | gawk '{print $1;}'`

#\rm $FAULTLIST
#foreach rup ( $RUPS )

#set REALS = `\ls $RUPDIR/$rup | gawk '{print $1;}'`
#foreach real ( $REALS )

#echo $RUPDIR/$rup/$real/${rup}_${real}.txt >> $FAULTLIST

#end
#end

foreach stat ( $STATS )

set SRC_CORDS = ../../data/SgtInfo/FdLocs/${SITE}.fdloc
#set SRC_CORDS = FdLocs/${stat}${SUFX}.fdloc
set OUTFILE = $OUTDIR/${SITE}.cordfile
#set OUTFILE = $OUTDIR/${stat}${SUFX}.cordfile

set XSRC = `gawk '{printf "%d\n",$1;}' $SRC_CORDS `
set YSRC = `gawk '{printf "%d\n",$2;}' $SRC_CORDS `

echo $XSRC $YSRC

gen_sgtgrid nx=$NX ny=$NY nz=$NZ h=$HH xsrc=$XSRC ysrc=$YSRC \
            ixmin=$IX_MIN ixmax=$IX_MAX iymin=$IY_MIN iymax=$IY_MAX \
	    izstart=$IZ_START izmax=$IZ_MAX \
	    radiusfile=$RADIUS_FILE outfile=$OUTFILE \
	    modellon=$MODEL_LON modellat=$MODEL_LAT modelrot=$MODEL_ROT \
	    faultlist=$FAULTLIST

end
